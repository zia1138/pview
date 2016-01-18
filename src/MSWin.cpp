#include <iomanip>
#include <fstream>
#include <math.h> 

#include "MSWin.hpp"
#include "xml.hpp"

#include <Qt>

// ================================================
// ================================================
// BEGIN - Drawing routines for main GUI. =========


// Draw peak data as pixels or as circles.
void LCMSView::draw_peaks(vector<lcms::peak *> &data, vector<lcms::ms2peak *> &data2, QPainter *painter) {
  double scaletime = double(width) / (t.max_time - t.min_time);
  double scalemz = double(height) / (d.max_mz - d.min_mz);
  double logImin = analysis->conf.logImin_disp, logImax = analysis->conf.logImax_disp;
  // Draw peaks with no MS/MS data.
  QPoint ip;
  for(size_t i = 0; i < data.size(); i++){
    lcms::peak *p = data[i];
    
    // Find color to use to show intensity.
    // NOTE: This code is repeated in XICView.
    double logI = log10(p->intensity);
    int level = 0;
    if(logI > logImax) level = 255; 
    else if( logI < logImin ) level = 0; 
    else level = (int)((logI - logImin) / (logImax - logImin) * 256.0);

    ip.rx() = (int)((p->retentionTime - t.min_time) * scaletime);
    ip.ry() = (int)((p->mz - d.min_mz) * scalemz );      
    // NOTE: Using setPixel is *much* faster than drawPoint.
    if(peak_size == 1) canvas->setPixel(ip, colors[level].rgb());
    else {
      QPen outline(colors[level]);
      QBrush fill(colors[level]);
      painter->setPen(outline); painter->setBrush(fill);
      painter->drawEllipse(ip.x() - peak_size, ip.y() - peak_size, 2 * peak_size, 2 * peak_size);
    }
  }

  // Draw MS/MS peaks on top in red.
  QColor red(255,0,0);
  QPen outline(red); QBrush fill(red); 
  painter->setPen(outline); 
  painter->setBrush(fill);
  for(size_t i = 0; i < data2.size(); i++){
    lcms::ms2peak *p = data2[i];
    ip.rx() = (int)((p->retentionTime - t.min_time) * scaletime);
    ip.ry() = (int)((p->mz - d.min_mz) * scalemz );      
    if(peak_size == 1) canvas->setPixel(ip, red.rgb());
    else painter->drawEllipse(ip.x() - peak_size, ip.y() - peak_size, 2 * peak_size, 2 * peak_size);
  }
}

// Display two types of query rectangles, useful for configurtion.
void LCMSView::draw_qrect(lcms::irun *irun, vector<lcms::peak *> &data, QPainter *painter) {
  if(display_qrect == FILTER_QRECT || display_qrect == NEIGHBOR_QRECT ) {

    double dmz;
    if(display_qrect == FILTER_QRECT) dmz = irun->conf.dmz; else dmz = irun->conf.peak_dmz;  

    // Draw a rectangle around each peak indicating the orthognal
    // range queries used to build the peak graph and filter noise.
    QPen outline(QColor(70,70,70)); 
    painter->setPen(outline); painter->setBrush(Qt::NoBrush);
    QRect rect;
    for(size_t i = 0; i < data.size(); i++){
      lcms::peak *p = data[i];
      float startTime, endTime; irun->get_range(startTime, endTime, p->retentionTime);
      int bx1, by1; to_canvas(startTime, p->mz - mass::ppm2da(p->mz, dmz) / 2, bx1, by1);
      int bx2, by2; to_canvas(endTime, p->mz + mass::ppm2da(p->mz, dmz) / 2, bx2, by2);
      rect.setX(bx1); rect.setY(by1);
      rect.setRight(bx2); rect.setBottom(by2);
      painter->drawRect(rect);
    }
  }
}

// Draws isotope labeled groups found in current view.
void LCMSView::draw_isotope_groups(lcms::irun *irun, QPainter *painter) {
  // Get isotope groups with current m/z and rt ranges.
  vector<lcms::isotope_group *> groups;
  irun->queryIsotopeGroups(groups, t.min_time, t.max_time, d.min_mz, d.max_mz);

  QPen outline(QColor(0,255,255)); 

  vector<lcms::xic*> xics;
  // Connects the ends of XICs to display isotope groups.
  for(size_t i = 0; i < groups.size(); i++) {
    // Set color based on log2 ratio vlaue.
    if(fabs(groups[i]->log2ratioHL) > 3) outline.setColor(QColor(255,0,0));
    else if(fabs(groups[i]->log2ratioHL) > 1) outline.setColor(QColor(0,255,255));
    else outline.setColor(QColor(128,128,128));
    painter->setPen(outline); 

    xics.clear();
    if(groups[i]->xics.size() > 0) {
      for(size_t x = 0; x < groups[i]->xics.size(); x++)
	if(groups[i]->xics[x] != NULL) xics.push_back(groups[i]->xics[x]);
    }
    else {
      if(groups[i]->light != NULL) xics.push_back(groups[i]->light);
      if(groups[i]->medium != NULL) xics.push_back(groups[i]->medium);
      if(groups[i]->heavy != NULL) xics.push_back(groups[i]->heavy);
    }

    for(size_t x = 1; x < xics.size(); x++) {
      double rt1L = xics[x-1]->start, rt2L = xics[x-1]->end, mzL = xics[x-1]->mz;
      double rt1H = xics[x]->start, rt2H = xics[x]->end, mzH = xics[x]->mz;
      int x1, y1; to_canvas(rt1L, mzL, x1, y1);
      int x2, y2; to_canvas(rt1H, mzH, x2, y2);
      painter->drawLine(x1,y1,x2,y2); 
      to_canvas(rt2L, mzL, x1, y1);
      to_canvas(rt2H, mzH, x2, y2);
      painter->drawLine(x1,y1,x2,y2); 
    }
  }
}

// Draw eXtracted Ion Chromatograms. 
void LCMSView::draw_XICs(lcms::irun *irun, QPainter *painter) {
  // Get XICs within current m/z range.
  vector<lcms::xic *> xics;
  irun->queryxics(xics, t.min_time, t.max_time, d.min_mz, d.max_mz);  

  // TODO: Change this so that  it also draws XICs that cut through the range!!! 
  QPen outline(QColor(128,0,128)); QBrush fill(QColor(128,0,128));
  painter->setBrush(fill);

  double scaletime = double(width) / (t.max_time - t.min_time);
  double scalemz = double(height) / (d.max_mz - d.min_mz);
  QColor red(255,0,255);

  for(unsigned int r = 0; r < xics.size(); r++){
    // Color differently if the XIC has charge data from an MS/MS ID.
    int ids = 0;
    if(xics[r]->ms2peaks.size() > 0) {
      ids = xics[r]->num_ids();
      if(ids > 0) {
	outline.setColor(QColor(0,255,255));
	fill.setColor(QColor(0,255,255));
      }
      else {
	if(xics[r]->type == lcms::xicMONOISO) {
	  outline.setColor(QColor(0,255,0));
	  fill.setColor(QColor(0,255,0));
	}
	else { outline.setColor(QColor(255,255,255)); fill.setColor(QColor(255,255,255));}
      }
    }
    else {
      if(xics[r]->type == lcms::xicMONOISO) { outline.setColor(QColor(0,100,0)); fill.setColor(QColor(0,100,0)); }
      else {
	if(xics[r]->type == lcms::xicISOTOPE) { outline.setColor(QColor(200,0,0)); fill.setColor(QColor(200,0,0)); }
	else { outline.setColor(QColor(128,128,128)); fill.setColor(QColor(128,128,128)); }
      }
    }

    painter->setPen(outline); painter->setBrush(fill);

    // Draw the actual XIC as a rectangle that scales with peak size. 
    int xL,yL; to_canvas(xics[r]->start, xics[r]->mz - mass::ppm2da(xics[r]->mz, analysis->conf.xic_width) / 2, xL,yL);
    int xR,yR; to_canvas(xics[r]->end, xics[r]->mz + mass::ppm2da(xics[r]->mz, analysis->conf.xic_width) / 2, xR, yR);

    if(ids > 0 && display_IDs) {
      vector<ms2::fragment *> ids;
      xics[r]->max_ms2()->max_ids(ids);
      
      QString info;
      for(size_t i = 0; i < ids.size(); i++) {
	info += analysis->get_meta(ids[i]).substr(0,20).c_str();
	if(i != ids.size() - 1)
	  info += " or ";
      }
      painter->drawText(xL,yL, info);
    }

    QRect rect(xL,yL,0,0); rect.setRight(xR); rect.setBottom(yR);
    painter->drawRect(rect);

    // Draw center of the XIC.
    outline.setColor(QColor(0,0,255)); painter->setPen(outline);
    fill.setColor(QColor(0,0,255)); painter->setBrush(fill);

    int x,y; to_canvas(xics[r]->retentionTime, xics[r]->mz, x, y);
    int dy = yR-yL;
    QRect rect_cent(x-dy,y-dy,0,0); rect_cent.setRight(x+dy); rect_cent.setBottom(y+dy);
    painter->drawEllipse(rect_cent);

    // Draws a line connecting all of the (rt, m/z) points contained
    // in the XIC.
    outline.setColor(QColor(255,255,255));
    painter->setPen(outline);
    for(unsigned int p = 1; p < xics[r]->chrom.size(); p++){
      double mz1 = xics[r]->chrom[p-1].mz;
      double retentionTime1 = xics[r]->chrom[p-1].retentionTime;
	  
      double mz2 = xics[r]->chrom[p].mz;
      double retentionTime2 = xics[r]->chrom[p].retentionTime;
      int x1, y1; to_canvas(retentionTime1, mz1, x1, y1);
      int x2, y2; to_canvas(retentionTime2, mz2, x2, y2);
      painter->drawLine(x1,y1,x2,y2);
    }

    outline.setColor(red); fill.setColor(red); 
    painter->setPen(outline); painter->setBrush(fill);
    QPoint ip;
    for(size_t k = 0; k < xics[r]->ms2peaks.size(); k++) {
      lcms::ms2peak &p = xics[r]->ms2peaks[k];
      ip.rx() = (int)((p.retentionTime - t.min_time) * scaletime);
      ip.ry() = (int)((p.mz - d.min_mz) * scalemz );      
      if(peak_size > 1) painter->drawEllipse(ip.x() - peak_size, ip.y() - peak_size, 2 * peak_size, 2 * peak_size);
    }
  }

  // Display query rectangles corresponding to grouping across
  // iruns/replicates/conditions.
  if(display_corr) {
    outline.setColor(QColor(0,200,200));
    painter->setPen(outline);

    for(unsigned int r = 0; r < xics.size(); r++){
      int x1 = 0, y1 = 0; 
      int x2 = 0, y2 = 0; 
      to_canvas(xics[r]->retentionTime - irun->conf.align_dtime / 2, 
		xics[r]->mz - mass::ppm2da(xics[r]->mz,irun->conf.xic_width) / 2, x1,y1);
      to_canvas(xics[r]->retentionTime + irun->conf.align_dtime / 2, 
		xics[r]->mz + mass::ppm2da(xics[r]->mz, irun->conf.xic_width) / 2, x2, y2);
      painter->drawLine(x1, y1, x2, y1);
      painter->drawLine(x2, y1, x2, y2);
      painter->drawLine(x2, y2, x1, y2);
      painter->drawLine(x1, y2, x1, y1);
    }    
  }

  if(display_group) {
    outline.setColor(QColor(0,200,200));
    painter->setPen(outline);

    for(unsigned int r = 0; r < xics.size(); r++){
      int x1 = 0, y1 = 0; 
      int x2 = 0, y2 = 0; 
      to_canvas(xics[r]->retentionTime - irun->conf.group_dtime / 2, 
		xics[r]->mz - mass::ppm2da(xics[r]->mz,irun->conf.xic_width) / 2, x1,y1);
      to_canvas(xics[r]->retentionTime + irun->conf.group_dtime / 2, 
		xics[r]->mz + mass::ppm2da(xics[r]->mz, irun->conf.xic_width) / 2, x2, y2);
      painter->drawLine(x1, y1, x2, y1);
      painter->drawLine(x2, y1, x2, y2);
      painter->drawLine(x2, y2, x1, y2);
      painter->drawLine(x1, y2, x1, y1);
    }    
  }


  if(display_isotopic) {
    painter->setBrush(Qt::NoBrush);
    vector<double> labels;

    /*for(size_t i = 0; i < irun->conf.isotopes.size(); i++) {
      if(irun->conf.isotopes[i].active)
	labels.push_back(irun->conf.isotopes[i].shift);
	}*/

    for(unsigned int r = 0; r < xics.size(); r++){
      if(xics[r]->charge > 0) {
	// Shows region where close-by XICs are filtered out.
	outline.setColor(QColor(0,200,0));
	painter->setPen(outline);
	int x1, y1; to_canvas(xics[r]->start, xics[r]->mz - mass::ppm2da(xics[r]->mz, irun->conf.xic_width), x1,y1);
	int x2, y2; to_canvas(xics[r]->end, xics[r]->mz + mass::ppm2da(xics[r]->mz,irun->conf.xic_width), x2, y2);
	painter->drawLine(x1, y1, x2, y1); painter->drawLine(x2, y1, x2, y2);
	painter->drawLine(x2, y2, x1, y2); painter->drawLine(x1, y2, x1, y1);

	// Shows query boxes used to find isotope groups.
	/*for(size_t k = 0; k < labels.size(); k++) {
	  outline.setColor(QColor(160,0,0));
	  painter->setPen(outline);
	  int charge = xics[r]->charge;
	  double isotopic_spacing = labels[k] / double(charge);
	  double mz_check = xics[r]->mz + isotopic_spacing;
 
	  int x1, y1; to_canvas(xics[r]->start, mz_check - mass::ppm2da(mz_check, irun->conf.isotolerence), x1,y1);
	  int x2, y2; to_canvas(xics[r]->end, mz_check + mass::ppm2da(mz_check,irun->conf.isotolerence), x2, y2);
           
	  painter->drawLine(x1, y1, x2, y1); painter->drawLine(x2, y1, x2, y2);
	  painter->drawLine(x2, y2, x1, y2); painter->drawLine(x1, y2, x1, y1);
	  }*/
      }
    }
  }

}

// Draws XIC groups.
void LCMSView::draw_grouped(QPainter *painter, vector<lcms::align_group *> &obs, QColor color) {
  QPen outline(color); QBrush fill(color);
  painter->setBrush(fill);
  painter->setPen(outline);    
  for(size_t i = 0; i < obs.size(); i++){
    int x,y; to_canvas(obs[i]->retentionTime, obs[i]->mz, x, y);

    for(size_t j = 0; j < obs[i]->xics.size(); j++){
      int xp,yp; to_canvas(obs[i]->xics[j]->retentionTime, obs[i]->xics[j]->mz, xp, yp);
      painter->drawLine(x,y,xp,yp);
      painter->drawEllipse(xp - 2 * peak_size, yp - 2 * peak_size, 4 * peak_size, 4 * peak_size);
    }
    painter->drawEllipse(x - 4 * peak_size, y - 4 * peak_size, 8 * peak_size, 8 * peak_size);
  }  
}

// Draw grouped observations.
void LCMSView::draw_grouped(QPainter *painter){
  // Get all observations within the current display window.
  vector<lcms::align_group *> obs; analysis->query(obs, t.min_time, t.max_time, d.min_mz, d.max_mz);
  draw_grouped(painter, obs, QColor(255,255,0));
}

// Draw graph connecting grouped/neighbor peaks.
void LCMSView::draw_graph(lcms::irun *irun, QPainter *painter) {
  // Get the edges in the query rectangle.
  QPen pen(QColor(255,255,0));
  painter->setPen(pen);
  vector< pair<lcms::peak*,lcms::peak*> > edges;
  irun->querypeakgraph(edges, t.min_time, t.max_time, d.min_mz, d.max_mz);
  for(unsigned int e = 0; e < edges.size(); e++){
    // Draw each edge.
    lcms::peak *p1 = edges[e].first;
    lcms::peak *p2 = edges[e].second;
    int x1,y1; to_canvas(p1->retentionTime, p1->mz, x1,y1);
    int x2,y2; to_canvas(p2->retentionTime, p2->mz, x2,y2);
    painter->drawLine(x1, y1, x2, y2);
  }
}

// Redraw entire LCMS data set using current zoom settings. TOOD: Put
// redraw in its own thread.
void LCMSView::redraw() {
  lcms::irun *irun = analysis->getIrun(condIndex, repIndex, irunIndex);

  // Use the QImage interface. 
  canvas->fill(QColor(0,0,0).rgb()); // black background.
  QPainter *painter = new QPainter(canvas);

  // Query peaks, raw or filtered.
  vector<lcms::peak *> data; 
  vector<lcms::ms2peak *> data2; 
  const double eps = 1e-5;
  if(display_raw) {
    irun->queryraw(data, t.min_time, t.max_time - eps, d.min_mz, d.max_mz - eps);
    irun->queryraw2(data2, t.min_time, t.max_time - eps, d.min_mz, d.max_mz - eps);
  }
  else {
    irun->queryfiltered(data, t.min_time, t.max_time - eps, d.min_mz, d.max_mz - eps);
    irun->queryfiltered2(data2, t.min_time, t.max_time - eps, d.min_mz, d.max_mz - eps);
  }

  // Draw query rectangles, if requested.
  draw_qrect(irun, data, painter);

  // Draw peak graph.
  if(display_graph && display_raw == false)  draw_graph(irun, painter);
  // Display XICs.
  if(display_XICs) {
    draw_XICs(irun, painter);
    draw_isotope_groups(irun, painter);
  }

  // Draw peak data.
  draw_peaks(data, data2, painter);

  // Draw symbols that show grouping of XICs across iruns/conditions.
  if(display_obs) draw_grouped(painter);

  delete painter;
}

void LCMSView::mouseReleaseEvent( QGraphicsSceneMouseEvent *) {
  if(selection != NULL) {
    scene()->removeItem(selection);
    QRectF rect = selection->rect().normalized();

    // Convert selection into m/z and retention time coordinates. 
    double rt1, mz1; from_canvas(rect.x(), rect.y(), rt1, mz1);
    double rt2, mz2; from_canvas(rect.right(), rect.bottom(), rt2, mz2);

    if(select_mode == ZOOM) { setmzrt(mz1, mz2, rt1, rt2); }
    else {
      lcms::irun *irun = analysis->getIrun(condIndex, repIndex, irunIndex);
      if(select_mode == MS2_SELECT) { 
	// Dispaly selected MS/MS peaks.
	vector<lcms::ms2peak *> peaks; 
	// Get the peaks in the selection region.
	if(display_raw) irun->queryraw2(peaks, rt1, rt2, mz1, mz2);
	else irun->queryfiltered2(peaks, rt1, rt2, mz1, mz2);
	// Fine the first one that has MS2 data.
	if(ms2win == NULL) ms2win = new MSMSWin(analysis);
	ms2win->show(); ms2win->raise(); ms2win->activateWindow();
	for(size_t i = 0; i < peaks.size(); i++)  ms2win->setPeak(peaks[i], irun->description);
      }
      else if(select_mode == XIC_SELECT) {
	// Display selected XICs.
	std::vector<lcms::xic *> xics;
	if(xicwin == NULL) xicwin = new XICWin(analysis->conf);
	xicwin->show(); xicwin->raise();  xicwin->activateWindow();
	irun->queryxics(xics, rt1, rt2, mz1, mz2);
	for(size_t i = 0; i < xics.size(); i++) xicwin->setXIC(xics[i], irun->description);
      }
      else if(select_mode == XIC_MS2_SELECT) {
	// Display MS/MS spectra grouped into an XIC.
	vector<lcms::xic *> xics; irun->queryxics(xics, rt1, rt2, mz1, mz2);
	if(ms2win == NULL) ms2win = new MSMSWin(analysis);
	ms2win->show();	ms2win->raise(); ms2win->activateWindow();
	for(size_t i = 0; i < xics.size(); i++) 
	  for(size_t j = 0; j < xics[i]->ms2peaks.size(); j++) 
	    ms2win->setPeak(&xics[i]->ms2peaks[j], irun->description);
      }
    }
    
    // Remove selection rectangle.
    delete selection;
    selection = NULL; select_mode = NO_SELECT;
  }
}

// END drawing routines ===========================
// ================================================
// ================================================

// Begin Main LC-MS data window.
MSWin::MSWin(QMainWindow *parent) : QMainWindow(parent) { 
  setupUi(this); 
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  graphicsView->setScene(&scene);

  // Initialize main data objects.
  view = NULL;
  analysis = NULL; 

  // Make column with protein group inforation larger.
  isotopeTree->header()->resizeSection(0, 550);
}

// Shortcut keyboard interface.
void MSWin::keyReleaseEvent(QKeyEvent *event) {
  if(view != NULL) {
    switch(event->key()) {
    case Qt::Key_A: view->cycleCorrRect(); break;
    case Qt::Key_P: view->toggleGroup(); break;

    case Qt::Key_Z: view->setZoomMode(); break;
    case Qt::Key_X: view->zoomOut(); break;
 
    case Qt::Key_G: view->toggleGraph(); break;
    case Qt::Key_F: view->toggleRawData(); break; 
    case Qt::Key_R: view->cycleRect(); break;

    case Qt::Key_1: view->setXICSelect(); break;
    case Qt::Key_2: view->setXICMS2Select(); break;
    case Qt::Key_3: view->setMS2Select(); break;

    case Qt::Key_S: view->toggleXICs(); break;
    case Qt::Key_D: view->toggleIDs(); break;
    case Qt::Key_I: view->toggleIsotopic(); break;
      
    case Qt::Key_Plus: view->incrementPeakSize(); break; 
    case Qt::Key_Minus: view->decrementPeakSize(); break; 

    case Qt::Key_O: view->toggleObs(); break;
    }
    event->accept();
  }
}

// Finds files with given suffix.
void MSWin::find_files(QString dir, vector<string> &base, vector<string> &files, QString suf) {
  // Iterate through directory.
  QDirIterator it(dir);
  while(it.hasNext()) {
    it.next();
    QFileInfo fileInfo = it.fileInfo();
    // Find any fasta files and save their full path.
    if(QString::compare(fileInfo.suffix(), suf, Qt::CaseInsensitive) == 0) {
      base.push_back(fileInfo.completeBaseName().toStdString());
      files.push_back(fileInfo.absoluteFilePath().toStdString());
    }
  }
}

// Collects mzXML files.
void collect_mzXML_files(QDir dir, vector<lcms::irun_file> &irunFiles) {
  QDirIterator it_dir(dir);
  while(it_dir.hasNext()) {
    it_dir.next();
    // Now only get files.
    QFileInfo info_file = it_dir.fileInfo(); if(info_file.isDir()) continue;

    // Save all mzXML files that were found.
    if(info_file.suffix() == "mzXML") {
      // We're now at DataSetXXX/ConditionY/RepZ/fileA.mzXML.
      lcms::irun_file irun;
      irun.description = info_file.baseName().toStdString();
      irun.filename = info_file.filePath().toStdString();
      irunFiles.push_back(irun);
    }
  }
  sort(irunFiles.begin(), irunFiles.end());
}

// Loads mzXML files and description information from subdirectories.
// Handles 3-level data structure.
void find_mzXML_files_3level(QString dir, lcms::analysis_files &analysis) {
  QDir dirInfo(dir);

  // Full analysis description extracted from directory name.
  analysis.description = dirInfo.dirName().toStdString();

  // Examine each subdirectory (these are the conditions).
  QDirIterator it_dir(dir);
  while(it_dir.hasNext()) { //  Level 1, the data set.
    // Iterate through the selected directory looking for subdirectories containing condition data.
    it_dir.next(); 
    QFileInfo info_subdir = it_dir.fileInfo();
    // Make sure it's a real directory.
    if(!info_subdir.isDir() || info_subdir.fileName() == ".." || info_subdir.fileName() == ".") continue;

    // We're now in DataSetXXXX/ConditionY.
    QDir subdir(info_subdir.filePath());
    lcms::condition_files condition; condition.description = subdir.dirName().toStdString();

    QDirIterator it_subdir(subdir);
    while(it_subdir.hasNext()) { // Level 2 a condition.
      it_subdir.next();       
      QFileInfo info_subsubdir = it_subdir.fileInfo();

      // Make sure it's a real directory.
      if(!info_subsubdir.isDir() || info_subsubdir.fileName() == ".." || info_subsubdir.fileName() == ".") continue;

      //  We're now in DataSetXXXX/ConditionY/RepZ.
      lcms::replicate_files replicate;
      QDir subsubdir(info_subsubdir.filePath());
      replicate.description = subsubdir.dirName().toStdString();

      // Look for reversed isotope labels flag.
      replicate.label_reversed = replicate.description.find("REVERSED") != string::npos;

      // For this replicate set, collect mzXML files
      // Level 3 a replicate set.
      collect_mzXML_files(subsubdir, replicate.irunFiles);

      // If there were files found for replicate add to condition.
      if(replicate.irunFiles.size() > 0) condition.replicateFiles.push_back(replicate);  
      
    }
    // If there were replicates found, add to analysis.
    if(condition.replicateFiles.size() > 0) {
      sort(condition.replicateFiles.begin(), condition.replicateFiles.end());
      analysis.conditionFiles.push_back(condition); 
    }
  }
  sort(analysis.conditionFiles.begin(), analysis.conditionFiles.end());
}

// Alternate single directory structure DataSet
void find_mzXML_files_1level(QString dir, lcms::analysis_files &analysis) {
  QDir dirInfo(dir); // Selected directory.
  // Full analysis description extracted from directory name.
  analysis.description = dirInfo.dirName().toStdString();

  lcms::condition_files conditionXXXX;
  conditionXXXX.description = analysis.description;
  // We have no replicate set subdirectory. Use XXX description.
  lcms::replicate_files replicateXXXX; replicateXXXX.description = "NA";
  replicateXXXX.label_reversed = false; // No reverse labeling allowed in two-level.

  collect_mzXML_files(dirInfo, replicateXXXX.irunFiles);

  if(replicateXXXX.irunFiles.size() > 0 ) {
    conditionXXXX.replicateFiles.push_back(replicateXXXX);
    analysis.conditionFiles.push_back(conditionXXXX);   
  }
}

// Alternate two-level directory structure. DataSet/ConditionX
void find_mzXML_files_2level(QString dir, lcms::analysis_files &analysis) {
  QDir dirInfo(dir); // Selected directory.
  // Full analysis description extracted from directory name.
  analysis.description = dirInfo.dirName().toStdString();
  // Examine each subdirectory (these are the conditions).
  QDirIterator it_dir(dir);
  while(it_dir.hasNext()) { // Level 1 data set.
    // Iterate through the selected directory looking for subdirectories containing condition data.
    it_dir.next(); 
    QFileInfo info_subdir = it_dir.fileInfo();
    // Make sure it's a real directory.
    if(!info_subdir.isDir() || info_subdir.fileName() == ".." || info_subdir.fileName() == ".") continue;

    // We're now in DataSetXXXX/ConditionY.
    QDir subdir(info_subdir.filePath());
    lcms::condition_files condition; condition.description = subdir.dirName().toStdString();

    // We have no replicate set subdirectory. Use XXX description.
    lcms::replicate_files replicateXXXX; replicateXXXX.description = "NA";
    replicateXXXX.label_reversed = false; // No reverse labeling allowed in two-level.

    // Level 2 condition subdirectory.
    // For this replicate set, collect mzXML files
    collect_mzXML_files(subdir, replicateXXXX.irunFiles);

    if(replicateXXXX.irunFiles.size() > 0) condition.replicateFiles.push_back(replicateXXXX);

    // If there was an "dummy" replicates set found w/ files, add to
    // analysis (no need to sort since there will only be one).
    if(condition.replicateFiles.size() > 0) analysis.conditionFiles.push_back(condition); 
  }
}

void MSWin::on_action_Recalibrate_triggered() {
  if(analysis == NULL) return;
  RecalWin recalwin(this); 
  if(recalwin.exec() == QDialog::Rejected) return;   
  
  lcms::recalmode_t recal;
  double winSz = 0;
  if(recalwin.recalMzClicked())   { recal = lcms::recalmodeMZ;   winSz = recalwin.getMzWin(); }
  if(recalwin.recalTimeClicked()) { recal = lcms::recalmodeTIME; winSz = recalwin.getTimeWin(); }

  // TODO: Make this multithreaded.
  for(int idx = 0; idx < analysis->num_iruns(); idx++) analysis->recalibrate_irun(recal, winSz, idx);
}

// Displays a dialog box that prompts the user for the data type used
// for conversion.
bool MSWin::profileDialog(bool &ms1_profile, bool &ms2_profile) {
  QMessageBox profileBox;
  profileBox.setIcon(QMessageBox::Question);
  profileBox.setWindowTitle("Centroid vs. Profile Data");
  profileBox.setText(tr("How was this data collected?"));
  profileBox.setInformativeText(tr("If the data was converted using Thermo2PVIEW, click all centroid. Otherwise, find out how these data were collected. PVIEW will not perform correctly without this information!"));
  QPushButton *centroidButton = profileBox.addButton(tr("All Centroid"), QMessageBox::ActionRole);
  QPushButton *ms1ms2Button   = profileBox.addButton(tr("MS1 Profile Only"), QMessageBox::ActionRole);
  QPushButton *profileButton  = profileBox.addButton(tr("All Profile"), QMessageBox::ActionRole);
  QPushButton *cancelButton   = profileBox.addButton(tr("Cancel"), QMessageBox::RejectRole);
  profileBox.exec();
  if(profileBox.clickedButton() == centroidButton) { ms1_profile = ms2_profile = false; } 
  else if(profileBox.clickedButton() == ms1ms2Button) { ms1_profile = true; ms2_profile = false; }
  else if(profileBox.clickedButton() == profileButton) { ms1_profile = ms2_profile = true; }
  else if(profileBox.clickedButton() == cancelButton) return false;
  return true;
}


class AnalysisThread : public QThread {
public:
  AnalysisThread(lcms::analysis *a_) { analysis = a_; }
protected:
  lcms::analysis *analysis;
  void run() { analysis->process(); }
};

// TODO: Move stuff in on_action_Open_triggered into here (minus the
// GUI crap).  This can then be used to drive PVIEW at the command
// line.
struct DriverThread : public QThread {
  QString dir;
  QString errorTitle, errorMessage;
  int errorStatus;
  DriverThread(QString dir) : dir(dir) { errorStatus = 0; }
protected:
  void run() {
    // Find mzXML files, separated by condition.  Map from condition name to a vector of mzXML files.
    lcms::analysis_files files; 
    // First try 3-level directory structure dataSet/ConditionX/RepSetX
    find_mzXML_files_3level(dir, files);
    // Now 2-level dataSet/ConditionX
    if(files.conditionFiles.size() == 0) find_mzXML_files_2level(dir, files);
    // And now 1-level structure. dataSet
    if(files.conditionFiles.size() == 0) find_mzXML_files_1level(dir, files);

    if(files.conditionFiles.size() == 0) {
      errorTitle = "Data Load Problem - No Files Found";
      errorMessage = "Did not find any mzXML files. Did you organize the files and directories correctly? Are you sure you selected the right directory?";
      return;
    }



  }
};

// Loads mzXML data, runs all processing algorithms.
// NOTE: All data processing starts here!!!!
// NOTE: All data processing starts here!!!!
void MSWin::on_action_Open_triggered() {
  QString dir = QFileDialog::getExistingDirectory(this, "Open Directory", "",
						  QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
  if(dir == "") return;

  // Find mzXML files, separated by condition.  Map from condition name to a vector of mzXML files.
  lcms::analysis_files files; 
  // First try 3-level directory structure dataSet/ConditionX/RepSetX
  find_mzXML_files_3level(dir, files);
  // Now 2-level dataSet/ConditionX
  if(files.conditionFiles.size() == 0) find_mzXML_files_2level(dir, files);
  // And now 1-level structure. dataSet
  if(files.conditionFiles.size() == 0) find_mzXML_files_1level(dir, files);

  // After all directory structures expended quit.
  if(files.conditionFiles.size() == 0) {
    QMessageBox::critical(this, tr("Data Load Problem - No Files Found"), 
     tr("Did not find any mzXML files. Did you organize the files and directories correctly? Are you sure you selected the right directory?"));
    return;
  }

  // Sort files and directories by their names.
  files.sort_all();   

  // Load any fasta files in directory. Use fasta files for DB search.
  vector<string> base_fasta, fasta; find_files(dir, base_fasta, fasta, "fasta");
  
  // Load any pepXML files if they exist.
  vector<string> base_pepXML, pepXML; find_files(dir + "/pepxml", base_pepXML, pepXML, "xml");

  if(fasta.size() == 0 && pepXML.size() == 0) {
    QMessageBox::warning(this, "Data Load", "Did not find any PepXML or FASTA files in the selected directory. Are you sure you selected the right directory?.");
  }

  bool ms2_profile = false, ms1_profile = false;
  if(profileDialog(ms1_profile, ms2_profile) == false) return;

  // Load pview.xml if it exists.
  QFile configXML(dir + "/pview.xml");
  if(configXML.exists()) xml::load_config(dir.toStdString() + "/pview.xml", conf);

  // Now, display configuration dialog. 
  MSConfigWin configure(conf, this);
  if(configure.exec() == QDialog::Rejected) return; // User hit cancel.

  // Save pview.xml configuration automatically.
  xml::save_config(dir.toStdString() + "/pview.xml", conf);

  // Update configuration parameters re profile mode.
  conf.ms2_profile = ms2_profile; conf.ms1_profile = ms1_profile;

  // Clear out previous data set.
  if(analysis != NULL) { delete analysis; analysis = NULL; }

  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));  

  QDateTime startTime = QDateTime::currentDateTime(); // Record start time of data processing.

  // Allocate slots for mzXML files. 
  analysis = new lcms::analysis(conf, files, base_fasta, fasta, pepXML);

  // Load and build protein fragment database, etc.
  loadDBs();

  // Load and process in each file (multi-threaded single irun processing).
  processIruns();  

  // Multi-threaded irun/irun alignment.
  alignIruns();

  // Process all combined replicates in one single thread.
  AnalysisThread *thread = new AnalysisThread(analysis);
  // Run analysis thread and display dialog.
  progressDialog("Processing replicates/conditions. Please wait.", (QThread *)thread);
  delete thread;
  QDateTime endTime = QDateTime::currentDateTime(); // Record end time of processing.

  QApplication::restoreOverrideCursor();

  // Display an information box with processing information.
  QString procTimeStr = "Run started at " + startTime.toString() + " and finished at " + endTime.toString() + ".";
  QMessageBox::information(this, tr("Processing time"), procTimeStr);

  // Set LCMSView initial size and position.
  if(view != NULL) { scene.removeItem(view); delete view; view = NULL; }

  // Use the frame width to adjust view's width and height.
  int fw = graphicsView->frameWidth();
  view = new LCMSView(analysis, graphicsView->width() - 2*fw - 1, graphicsView->height() - 2*fw - 1, mzLabel, rtLabel);
  // If algorithms where not run, then toggle display to show all raw data.
  if(conf.run_algorithms == false) view->toggleRawData();
  scene.addItem(view);
  graphicsView->setSceneRect(0,0,graphicsView->width() - 2*fw-1, graphicsView->height() - 2*fw-1);

  // Empty lists on side of GUI.
  listCond->clear(); listRep->clear(); listIrun->clear();

  // Populate the conditions and replicates list.  Make sure 0th
  // replicate and 0th condition are selected!!!  NOTE: We only need
  // to poulate conditions because the select event populates the
  // other lists.
  for(int condIndex = 0; condIndex < analysis->numConds(); condIndex++) {
    QString condition = analysis->getCond(condIndex)->description.c_str();
    listCond->addItem(condition);
  }

  // Erase "memory" of selected replicate across conditions.
  selectedRepIndex.resize(analysis->numConds());
  // Reset the "memory" to the first replicate.
  for(size_t i = 0; i < selectedRepIndex.size(); i++) selectedRepIndex[i] = 0;

  // Erase "memory" of selected irun across conditions.
  selectedIrunIndex.resize(analysis->numConds());
  for(int condIndex = 0; condIndex < analysis->numConds();  condIndex++) {
    selectedIrunIndex[condIndex].resize(analysis->numReps(condIndex));
    // Set to 0th irun by default.
    for(size_t i = 0; i < selectedIrunIndex[condIndex].size(); i++)
      selectedIrunIndex[condIndex][i] = 0;
  }

  // Re-initialize positions in both lists.
  listCond->setCurrentRow(0); 

  // If in isotope mode, populate the isotope tree.
  isotopeTree->clear(); // Clear tree items.
  xicTree->clear();
  item2group.clear(); // Also clear mapping from item to isotope group.

  // Populate tree view with protein groups.
  if(conf.quantmode == quantISOPAIR || conf.quantmode == quantISOTRIPLE || conf.quantmode == quantAcetylSTOICH) populateIsotopeTree();
  else if(conf.quantmode == quantALIGN) populateXICTreeAlign();
  else if(conf.quantmode == quantXICBASED) populateXICTree();

}

QTreeWidgetItem * MSWin::accessionSubtree(vector<int> &idxs) {
  QTreeWidgetItem *acc = new QTreeWidgetItem(QStringList() << "Proteins");
  for(size_t j = 0; j < idxs.size(); j++) {
    QString meta = analysis->get_meta(idxs[j]).c_str();
    // Add new-lines to meta.
    int i = 0; while(i < meta.length()) { if(i != 0 && i % 50 == 0) { meta.insert(i, '\n'); }  i++; }
    acc->addChild(new QTreeWidgetItem(QStringList(meta)));
  }
  return acc;
}

struct xic_index { 
  bool operator () (lcms::xic* A, lcms::xic* B) { 
    return A->condIndex < B->condIndex && A->repIndex < B->repIndex && A->irunIndex < B->irunIndex; 
  } 
};

QTreeWidgetItem *MSWin::xicSubtree(string seq, int charge, vector<lcms::xic *> &xics) {
  // Use sequence and charge information to list XICs.
  QString fragStr = QString(seq.c_str()) + ", +" + QString::number(charge);
  QTreeWidgetItem *xicroot = new QTreeWidgetItem(QStringList() << fragStr);

  // Sort XICs based on condIndex, repIndex, irunIndex.
  sort(xics.begin(), xics.end(), xic_index());

  for(size_t i = 0; i < xics.size(); i++) {
    ms2::ms2peak *p = xics[i]->max_ms2(); if(p == NULL) continue;
    string seq = analysis->get_seqS(p);

    QString xicStr = QString(seq.c_str()) + "\n";

    xicStr +=  
      QString::number(xics[i]->quant, 'f', 2) + ":(" + QString::number(xics[i]->retentionTime) + "," + 
      QString::number(xics[i]->mz) + ")" + "\n";

    int condIndex = xics[i]->condIndex;
    int repIndex = xics[i]->repIndex;
    int irunIndex = xics[i]->irunIndex;

    QString condStr = analysis->getCond(condIndex)->description.c_str();
    QString repStr = analysis->getRep(condIndex,repIndex)->description.c_str();
    QString irunStr = analysis->getIrun(condIndex,repIndex,irunIndex)->description.c_str();
    
    QString sourceStr = condStr + "/" + repStr + "/" + irunStr;

    QTreeWidgetItem *xicItem = new QTreeWidgetItem(QStringList() << (xicStr + ":" + sourceStr));

    // Allow double click to view XIC.
    item2xic[xicItem] = xics[i];

    xicroot->addChild(xicItem);
  }
  return xicroot;
}


// Used to sort groups based on condition information.
struct group_index { 
  bool operator () (lcms::isotope_group* A, lcms::isotope_group* B) { 
    return A->condIndex < B->condIndex && A->repIndex < B->repIndex && A->irunIndex < B->irunIndex; 
  } 
};

void MSWin::populateIsotopeGroups(QTreeWidgetItem *frags, vector<lcms::isop_seqgroup> &seq_groups) {
  for(size_t s = 0; s < seq_groups.size(); s++) {
    QTreeWidgetItem *seq = new QTreeWidgetItem(QStringList() << QString(seq_groups[s].seq.c_str()));
    frags->addChild(seq);
    vector<lcms::isotope_group *> &groups = seq_groups[s].groups;
    sort(groups.begin(), groups.end(), group_index());
    for(size_t j = 0; j < groups.size(); j++) {

      // Extract condition + replicate + irun number
      int condIndex = groups[j]->condIndex;
      int repIndex = groups[j]->repIndex;
      int irunIndex = groups[j]->irunIndex;

      QString condStr = analysis->getCond(condIndex)->description.c_str();
      QString repStr = analysis->getRep(condIndex,repIndex)->description.c_str();
      QString irunStr = analysis->getIrun(condIndex,repIndex,irunIndex)->description.c_str();
    
      QString sourceStr = condStr + "/" + repStr + "/" + irunStr;

      lcms::xic *x = NULL;
      if(groups[j]->heavy != NULL) x = groups[j]->heavy;
      if(groups[j]->medium != NULL) x = groups[j]->medium;
      if(groups[j]->light != NULL) x = groups[j]->light; // We want to display the light XIC.

      if(groups[j]->xics.size() > 0) {
	for(size_t ii = 0; ii < groups[j]->xics.size(); ii++)  {
	  if(groups[j]->xics[ii] != NULL)
	    x = groups[j]->xics[ii];
	}
      }

      QString groupStr;
      if(x != NULL) 
	groupStr = "(" + QString::number(x->retentionTime) + "," + 
	  QString::number(x->mz) + "):" + sourceStr + ", +" + QString::number(x->charge);

      QString log2ratioHL = "NA", log2ratioML = "NA", log2ratioHM = "NA";
      if(groups[j]->heavy && groups[j]->light)  log2ratioHL = QString::number(groups[j]->log2ratioHL);
      if(groups[j]->medium && groups[j]->light) log2ratioML = QString::number(groups[j]->log2ratioML);
      if(groups[j]->heavy && groups[j]->medium) log2ratioHM = QString::number(groups[j]->log2ratioHM);
      
      if(log2ratioHL == "NA" && log2ratioML == "NA" && log2ratioHM == "NA") {
	if(groups[j]->heavy)  { log2ratioHL = "HEAVY"; log2ratioHM = "HEAVY"; }
	if(groups[j]->medium) { log2ratioML = "MEDIUM"; log2ratioHM = "MEDIUM"; }
	if(groups[j]->light)  { log2ratioHL = "LIGHT"; log2ratioHL = "LIGHT"; }
      }

      // Add item with column 1 with group info and column 2 with log2 ratio.
      QTreeWidgetItem *ip = 
	new QTreeWidgetItem(QStringList() << groupStr << log2ratioHL << log2ratioML << log2ratioHM);

      seq->addChild(ip);
      item2group[ip] = groups[j]; // Save in map to allow double click navigation.
    }
  }
}

void MSWin::populateIsotopeTree() {
  QTreeWidgetItem *root = new QTreeWidgetItem(QStringList("Protein Groups"));
  isotopeTree->addTopLevelItem(root);
  for(int groupIdx = 0; groupIdx < analysis->numProteinGroups(); groupIdx++) {
    lcms::isop_pgroup *group = analysis->isotopeGroupProteinGroup(groupIdx);
    QString groupStr = "Protein Group " + 
      QString::number(groupIdx + 1) + ", N = " + QString::number(group->groups.size());
    QTreeWidgetItem *item = new QTreeWidgetItem(QStringList() << groupStr);

    QTreeWidgetItem *acc = accessionSubtree(group->acc);

    QTreeWidgetItem *frags = new QTreeWidgetItem(QStringList("Fragments"));
    populateIsotopeGroups(frags, group->seq_groups);

    QTreeWidgetItem *shared_frags = NULL;

    if(group->seq_groups_shared.size() > 0) {
      shared_frags = new QTreeWidgetItem(QStringList("Shared Fragments"));
      populateIsotopeGroups(shared_frags, group->seq_groups_shared);
    }

    item->addChild(acc);
    item->addChild(frags);
    if(shared_frags != NULL) item->addChild(shared_frags);
    root->addChild(item);
  }
}

void MSWin::populateXICTreeAlign() {
  QTreeWidgetItem *root = new QTreeWidgetItem(QStringList("Protein Groups"));
  xicTree->addTopLevelItem(root);

  for(int groupIdx = 0; groupIdx < analysis->numProteinGroups(); groupIdx++) {
    QString groupStr = "Protein Group " + QString::number(groupIdx + 1);
    QTreeWidgetItem *item = new QTreeWidgetItem(QStringList() << groupStr); 

    lcms::align_pgroup *proteinGroup = analysis->alignProteinGroup(groupIdx);

    QTreeWidgetItem *acc = accessionSubtree(proteinGroup->acc);
    item->addChild(acc);  // Proteins -> Group i -> Accessions

    QTreeWidgetItem *groups = new QTreeWidgetItem(QStringList() << "Fragments");
    // Note: AlignGroups should be sorted on sizea nd score.
    vector<lcms::align_group *> &alignGroups = proteinGroup->align_groups;

    for(size_t i = 0; i < alignGroups.size(); i++) {
      lcms::align_group *alignGroup = alignGroups[i];
      lcms::ms2peak *p = alignGroup->max_ms2();

      if(p == NULL) continue;

      // Each align group has a sequence and a charge read off of the
      // maximum scoring MS/MS spectrum.
      string seq = analysis->get_seqS(p); int charge = p->charge;

      QTreeWidgetItem *xic = xicSubtree(seq, charge, alignGroup->xics);
      groups->addChild(xic);
    }
    item->addChild(groups); // Proteins -> Group i -> Groups

    QTreeWidgetItem *groups2 = new QTreeWidgetItem(QStringList() << "Shared Fragments");
    vector<lcms::align_group *> &sharedAlignGroups = proteinGroup->shared_align_groups;
    for(size_t i = 0; i < sharedAlignGroups.size(); i++) {
      lcms::align_group *sharedAlignGroup = sharedAlignGroups[i];
      lcms::ms2peak *p = sharedAlignGroup->max_ms2();

      if(p == NULL) continue;

      // Each align group has a sequence and a charge read off of the
      // maximum scoring MS/MS spectrum.
      string seq = analysis->get_seqS(p); int charge = p->charge;

      QTreeWidgetItem *xic = xicSubtree(seq, charge, sharedAlignGroup->xics);
      groups2->addChild(xic);
    }
    // List shared fragments.
    if(sharedAlignGroups.size() > 0)
      item->addChild(groups2); // Proteins -> Group i -> Groups
    else
      delete groups2;

    root->addChild(item); // Proteins -> Group i
  }
  
}

void MSWin::populateXICTree() {
  QTreeWidgetItem *root = new QTreeWidgetItem(QStringList("Protein Groups"));
  xicTree->addTopLevelItem(root);  
  for(int groupIdx = 0; groupIdx < analysis->numProteinGroups(); groupIdx++) {
    QString groupStr = "Protein Group " + QString::number(groupIdx + 1);
    QTreeWidgetItem *group = new QTreeWidgetItem(QStringList() << groupStr);

    lcms::xic_pgroup *proteinGroup = analysis->xicProteinGroup(groupIdx);
    QTreeWidgetItem *acc = accessionSubtree(proteinGroup->acc);
   
    QTreeWidgetItem *fragments = new QTreeWidgetItem(QStringList() << "Fragments");

    vector<lcms::seq_group> &seqGroups = proteinGroup->seq_groups;
    for(size_t i = 0; i < seqGroups.size(); i++) {
      QTreeWidgetItem *xic = xicSubtree(seqGroups[i].seq, seqGroups[i].charge, seqGroups[i].xics);
      fragments->addChild(xic);
    }

    QTreeWidgetItem *fragments_shared = new QTreeWidgetItem(QStringList() << "Shared Fragments");

    vector<lcms::seq_group> &seqGroupsShared = proteinGroup->seq_groups_shared;
    for(size_t i = 0; i < seqGroupsShared.size(); i++) {
      QTreeWidgetItem *xic = xicSubtree(seqGroupsShared[i].seq, seqGroupsShared[i].charge, seqGroupsShared[i].xics);
      fragments_shared->addChild(xic);
    }
    
    group->addChild(acc);  // Proteins -> Group i -> Accessions
    group->addChild(fragments); // Proteins -> Group i -> Groups

    if(seqGroupsShared.size() > 0) group->addChild(fragments_shared);
    else delete fragments_shared;

    root->addChild(group); // Proteins -> Group i
  }
}


// Thread for processing indivdual iruns.
class IrunThread : public QThread {
public:
  IrunThread(lcms::analysis *a_) { analysis = a_; }
  void addIrun(int idx) { irun_idxs.push_back(idx); }
protected:
  vector<int> irun_idxs; // Iruns assigned to thread.
  lcms::analysis *analysis;
  void run() { for(size_t i = 0; i < irun_idxs.size(); i++) analysis->process_irun(irun_idxs[i]); }
};

// Launches separate threads to process each replicate set.
void MSWin::processIruns() {
  // Allocate threads.
  vector<IrunThread *> threads(conf.threads);
  for(int i = 0; i < conf.threads; i++) threads[i] = new IrunThread(analysis);

  // Divide the replicates evenly among the threads.
  int T = 0;
  for(int irunIndex = 0; irunIndex < analysis->num_iruns(); irunIndex++) {
    // Give file lists and replicate sets to each thread.
    threads[T]->addIrun(irunIndex);
    T++; if(T >= conf.threads) T = 0;
  }

  // Run threads, display process dialog.
  vector<QThread *> threads2;
  for(size_t i = 0; i < threads.size(); i++) threads2.push_back(threads[i]);
  progressDialog("Loading and processing data. Please wait.", threads2);

  // Free up memory used by threads.
  for(int i = 0; i < conf.threads; i++) delete threads[i];
}


class AlignThread : public QThread {
public:
  AlignThread(lcms::analysis *a_, int refIdx_) { analysis = a_; refIdx = refIdx_; }
  void addIrun(int idx) { irun_idxs.push_back(idx); }
protected:
  int refIdx;  vector<int> irun_idxs; // Iruns assigned to thread.
  lcms::analysis *analysis;
  void run() { for(size_t i = 0; i < irun_idxs.size(); i++) analysis->align_irun(irun_idxs[i], refIdx); }
};


void MSWin::alignIruns() {
  bool ok_modes = conf.quantmode == quantALIGN;
  if(ok_modes == false) return;
  // Pick reference irun.
  int refIdx = analysis->selectAlignReference();
  cout << "align reference is " << analysis->iruns.at(refIdx)->description << endl;

  vector<AlignThread *> threads(conf.threads);
  for(int i = 0; i < conf.threads; i++) threads[i] = new AlignThread(analysis, refIdx);

  // Divide the replicates evenly among the threads.
  int T = 0;
  for(int irunIndex = 0; irunIndex < analysis->num_iruns(); irunIndex++) {
    if(irunIndex == refIdx) continue;
    // Give file lists and replicate sets to each thread.
    threads[T]->addIrun(irunIndex);
    T++; if(T >= conf.threads) T = 0;
  }
  // Run threads, display process dialog.
  vector<QThread *> threads2;
  for(size_t i = 0; i < threads.size(); i++) threads2.push_back(threads[i]);
  progressDialog("Aligning instrument runs. Please wait.", threads2);

  // Free up memory used by threads.
  for(int i = 0; i < conf.threads; i++) delete threads[i];
}

void MSWin::progressDialog(QString message, vector<QThread *> threads) {
  // Start all of the threads.
  for(int i = 0; i < (int)threads.size(); i++) threads[i]->start();
  /*QProgressDialog progress(message, 0,0, 0, this);
  progress.setWindowTitle("Pton LC-MS/MS Data Analyzer");
  progress.setModal(true); progress.show();*/
  QApplication::sendPostedEvents();
  while(true) {
    // Wait until both loader threads complete.
    bool done = true;
    for(size_t i = 0; i < threads.size(); i++) { done = done && threads[i]->wait(100); QApplication::sendPostedEvents(); }
    if(done) break;
  }
}

// Thread loads fasta file with AA sequence database.
// Thread also loads any AMT databases.
class AnalysisLoadThread : public QThread {
public:
  AnalysisLoadThread(lcms::analysis *a_) { analysis = a_; }
protected:
  lcms::analysis *analysis;
  void run() { analysis->load_protein_db(); }
};

// Loads and builds protein fragment database.
void MSWin::loadDBs() {
  AnalysisLoadThread *thread = new AnalysisLoadThread(analysis);
  progressDialog("Building forward and reverse fragment databases.", (QThread*)thread);
  delete thread;
}

// Output quantification values in TXT format.  TODO: Output this data
// in a separate thread so GUI doesn't hang.
void MSWin::on_actionSave_SILAC_TXT_triggered() {
  if(analysis == NULL) return;

  QString fileName = QFileDialog::getSaveFileName(this, tr("Save Tab Delimiated"), "", tr("TXT (*.txt)"));
  if(conf.quantmode == quantAcetylSTOICH) {
    QFileInfo file_info(fileName); 
    QString dir = file_info.dir().absolutePath();

    // Save raw isotope pairs with protein mappings, no protein groups.
    QString isodata_file = dir + "/" + file_info.baseName() + "_stoich.txt";
    analysis->save_stoich_data(isodata_file.toStdString());
    return;
  }  
  if(!(conf.quantmode == quantISOPAIR || conf.quantmode == quantISOTRIPLE)) return;
  if(fileName != "") {
    QFileInfo file_info(fileName); 
    QString dir = file_info.dir().absolutePath();
    // Save all isotope groups found, grouped by protein.
    analysis->save_isotope_groups(fileName.toStdString());

    // Save raw isotope pairs with protein mappings, no protein groups.
    QString isodata_file = dir + "/" + file_info.baseName() + "_isodata.txt";
    analysis->save_isotope_data(isodata_file.toStdString());


    // Save per-protein internal correlations. Used for group1 vs group2 plots.
    QString internal_file = dir + "/" + file_info.baseName() + "_internal.txt";
    analysis->save_isotope_groups_internal(internal_file.toStdString());

    // Save per-peptide internal correlations. Use for group1 vs group2 for peptides seen several times.
    QString internal_file2 = dir + "/" + file_info.baseName() + "_internal_pep.txt";
    analysis->save_isotope_groups_internal_pep(internal_file2.toStdString());

    // Save correlations between replicates. Used for replicate 1 vs replicate 2 correlations.
    QString corr_file = dir + "/" + file_info.baseName() + "_corr.txt";
    analysis->save_isotope_groups_corr(corr_file.toStdString());

    // Save mass errors.
    QString mass_error_file = dir + "/" + file_info.baseName() + "_mass_error_ms1.txt";
    analysis->saveMassErrorsMS1(mass_error_file.toStdString());

    QString mass_error2_file = dir + "/" + file_info.baseName() + "_mass_error_ms2.txt";
    analysis->saveMassErrorsMS2(mass_error2_file.toStdString());

    // Save ratios in per-protein tabular output.
    QString table_file = dir + "/" + file_info.baseName() + "_table.txt";
    analysis->save_isotope_groups_table(table_file.toStdString());

    // Save ratios in per-peptide tabular output
    QString table_file2 = dir + "/" + file_info.baseName() + "_table_pep.txt";
    analysis->save_isotope_groups_table_pep(table_file2.toStdString());
  }
}

// Saves full analysis to a TXT file.
void MSWin::on_actionSave_TXT_triggered() {
  if(analysis == NULL) return;
  if(! (conf.quantmode == quantXICBASED || conf.quantmode == quantALIGN) ) return;
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save Tab Delimiated"), "", tr("TXT (*.txt)"));
  if(fileName != "") {
    QFileInfo file_info(fileName); 
    QString dir = file_info.dir().absolutePath();
    QString table_peptides = dir + "/" + file_info.baseName() + "_peptides.txt";
    QString internal_file = dir + "/" + file_info.baseName() + "_internal.txt";
    QString mass_error_file = dir + "/" +file_info.baseName() + "_mass_error.txt";
    QString mass_error_file_ms2 = dir + "/" + file_info.baseName() + "_mass_error_ms2.txt";
    QString allid_file = dir + "/" + file_info.baseName() + "_allids.txt";
    if(conf.quantmode == quantALIGN) {
      // Save table of Aligned XICs by protein.
      analysis->saveAlignTable(fileName.toStdString());
      // Save table of Aligned XICs by peptide.
      analysis->saveAlignTablePeptides(table_peptides.toStdString());
      // Save table with internal correlations between all n * (n-1) / 2 groups of conditions.
      analysis->saveAlignTableInternal(internal_file.toStdString());
      // Save precursor mass errors.
      analysis->saveMassErrorsMS1(mass_error_file.toStdString());
      // Save MS/MS product ion mass errors.
      analysis->saveMassErrorsMS2(mass_error_file_ms2.toStdString());      
    }
    else if(conf.quantmode == quantXICBASED) {
      analysis->saveXICBasedTable(fileName.toStdString());
      analysis->saveXICBasedTablePeptides(table_peptides.toStdString());
      analysis->saveXICBasedTableInternal(internal_file.toStdString());
      analysis->saveMassErrorsMS1(mass_error_file.toStdString());      
      analysis->saveMassErrorsMS2(mass_error_file_ms2.toStdString());      
      analysis->saveXIC_allids(allid_file.toStdString());            
    }
  }
}

void MSWin::on_actionSave_Align_XICs_TXT_triggered() {
  if(analysis == NULL) return;
  if(conf.quantmode != quantALIGN) return;
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save Tab Delimiated"), "", tr("TXT (*.txt)"));
  if(fileName != "") analysis->saveLabelFreeMS1(fileName.toStdString());
}


void MSWin::on_isotopeTree_itemDoubleClicked(QTreeWidgetItem * item, int /*column*/ ) {
  if(item2group.find(item) == item2group.end()) return;

  lcms::isotope_group *group = item2group[item];

  vector<lcms::xic *> xics;
  if(group->light != NULL) xics.push_back(group->light);
  if(group->medium != NULL) xics.push_back(group->medium);
  if(group->heavy != NULL) xics.push_back(group->heavy);

  for(size_t i = 0; i < group->xics.size(); i++) {
    if(group->xics[i] != NULL) xics.push_back(group->xics[i]);
  }
  // Update the list boxes.
  if(xics.size() >= 2) {
    lcms::xic *light = xics[0], *heavy = xics[xics.size() - 1]; 
    listCond->setCurrentRow(light->condIndex);
    listRep->setCurrentRow(light->repIndex);
    listIrun->setCurrentRow(light->irunIndex);
    view->setmzrt(light->mz - 5, heavy->mz + 5,
		  light->retentionTime - 500, light->retentionTime + 500);
  }
  if(xics.size() == 1) {
    lcms::xic *cur = xics[0];
    listCond->setCurrentRow(cur->condIndex);
    listRep->setCurrentRow(cur->repIndex);
    listIrun->setCurrentRow(cur->irunIndex);
    view->setmzrt(cur->mz - 5, cur->mz + 5,
		  cur->retentionTime - 500, cur->retentionTime + 500);

  }

  tabMainDisplay->setCurrentIndex(0);
}


void MSWin::on_xicTree_itemDoubleClicked(QTreeWidgetItem * item, int /*column*/ ) {
  if(item2xic.find(item) == item2xic.end()) return;

  lcms::xic *xic = item2xic[item];

  // Update the list boxes.
  listCond->setCurrentRow(xic->condIndex); listRep->setCurrentRow(xic->repIndex); listIrun->setCurrentRow(xic->irunIndex);

  // Update the current zoom region.
  view->setmzrt(xic->mz - 5, xic->mz + 5,
		xic->retentionTime - 500, xic->retentionTime + 500);

  tabMainDisplay->setCurrentIndex(0);
}

// If user changed the condition row, update condition. 
void MSWin::on_listCond_currentRowChanged(int newCondIndex) {
  if(newCondIndex == -1 || view == NULL) return;
  listRep->clear();

  // Update the replicate list.
  for(int repIndex = 0; repIndex < analysis->numReps(newCondIndex); repIndex++) {
    QString replicate = analysis->getRep(newCondIndex, repIndex)->description.c_str();
    listRep->addItem(replicate);
  }

  // Update the condition number, do not update the view/redraw.
  view->setCond(newCondIndex); // Set new condition number.

  // Change the replicate, and update the view.
  listRep->setCurrentRow(selectedRepIndex[newCondIndex]);  
}

// If user changed the replicate row, update the replicate.  This code
// is arranged so it works correclty set setCurrentRow() call in
// listCond_currentRowChanged().
void MSWin::on_listRep_currentRowChanged(int newRepIndex) {
  if(newRepIndex == -1 || view == NULL) return;

  listIrun->clear();
  // Update irun list.
  for(int irunIndex = 0; irunIndex < analysis->numIruns(view->getCond(), newRepIndex); irunIndex++) {
    QString irun = analysis->getIrun(view->getCond(), newRepIndex, irunIndex)->description.c_str();
    listIrun->addItem(irun);
  }

  // Update replicate/fraction number, update the display in the view.
  view->setRep(newRepIndex);
  selectedRepIndex[view->getCond()] = newRepIndex; // Save memory of current replicate number.

  // Change the irun back to "memory" value.
  listIrun->setCurrentRow(selectedIrunIndex[view->getCond()][newRepIndex]);
}

// User changed the selected irun.
void MSWin::on_listIrun_currentRowChanged(int newIrunIndex) {
  if(newIrunIndex == -1 || view == NULL) return;
  // Save memory of selected irun.
  selectedIrunIndex[view->getCond()][view->getRep()] = newIrunIndex;
  // Update irun number in view.
  view->setIrun(newIrunIndex);
  // Redraw and update!
  view->redraw(); view->update(); 
}

// Application start point.
int main(int argc, char *argv[]) {
  // Launch main application.
  if(argc == 1) {
    QApplication app(argc, argv);
    MSWin pview; pview.show();
    return app.exec();
  }
  else {
    vector< vector<int> > res;
    /*util::multi_comb mc(8, 3, res);
      }*/
    vector<int> L;
    L.push_back(2);
    L.push_back(2);
    L.push_back(2);
    util::cross(L, res);
    for(size_t r = 0; r < res.size(); r++) {
      for(size_t i = 0; i < res[r].size(); i++) cout << res[r][i] << " "; cout << endl;
    }
    cout << "running command line" << endl;
  }
}


