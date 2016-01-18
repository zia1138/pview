#include <iostream>
#include "MSMSWin.hpp"
#include <Qt>

using namespace std;


// Begin - GUI for selected MS/MS fragmentation spectra.
// Draw MS/MS spectrum.
void MSMSView::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *) {
  QColor black(0,0,0);
  // Use black background for contrast.
  painter->setPen(QPen(black));
  painter->setBrush(QBrush(black));
  painter->drawRect(0,0, width, height);

  if(peak == NULL) return;  
  if(peak->ms2.size() == 0) return;

  QColor green(0,255,0), blue(0,0,255), magenta(255,0,255), white(255,255,255);
  QColor red(255,0,0), gray(150,150,150), darkmagenta(150,0,150);

  // Display MS/MS spectrum.
  double smz = width / (maxmz - minmz);
  double sI = height / maxI;

  // Display assigned b-ion y-ion ladder.
  if(fragment != NULL && ms2db != NULL) {
    vector<double> ladder;
    string seq;
    analysis->fasta.get_seq_nomod(fragment, seq);

    // The min_mz2 range for the ladder can be different for different mzXML files.
    ms2db->build_ladder(peak->activation, ladder, seq, *mod, peak->charge, peak->neutral_mass(), min_mz2, max_mz2);

    // Sort ladder from smallest m/z to largest m/z.
    sort(ladder.begin(), ladder.end());

    // Note that peak->ms2 is assumed to be sorted decreasing by
    // intensity.
    vector<bool> ladder_matched; vector<int> peak_matched;
    ms2db->match(peak->ms2, peak_matched, ladder, ladder_matched);

    // Draw theoretical peaks. Switch colors if theoretical ladder
    // peak is matched or not.
    for(size_t i = 0; i < ladder.size(); i++) {
      if(ladder_matched[i]) painter->setPen(QPen(green));
      else painter->setPen(QPen(blue));
      double x = (ladder[i] - minmz) * smz;
      painter->drawLine(QPointF(x, height), QPointF(x, 0));
    }

    for(size_t m = 0; m < mod->sites.size(); m++) {
      ladder.clear();
      ms2db->build_ladder(peak->activation, ladder, seq, *mod, peak->charge, peak->neutral_mass(), 
			  min_mz2, max_mz2, mod->sites[m].pos);
      painter->setPen(QPen(red));
      for(size_t i = 0; i < ladder.size(); i++) {
	double x = (ladder[i] - minmz) * smz;
	painter->drawLine(QPointF(x, height/2), QPointF(x, 0));
      }
    }

    // Draw observed data.
    vector<ms2::peak2> &ms2 = peak->ms2;
    for(size_t i = 0; i < ms2.size(); i++) {
      lcms::peak2 &p = ms2[i];
      if(p.mz2 < maxmz && p.mz2 > minmz) {
 	double x = (p.mz2 - minmz) * smz;
	// Use magenta for top-k peaks.
	if((int)i < 5) {
	  // Use darker color for unmatched peaks.
	  if(peak_matched[i] >= 0) painter->setPen(QPen(magenta));
	  else painter->setPen(QPen(darkmagenta));
	}
	else {
	  // Use darker color for unmatched peaks.
	  if(peak_matched[i] >= 0) painter->setPen(QPen(white));
	  else painter->setPen(QPen(gray));
	}
	// Draw peak.
 	painter->drawLine(QPointF(x, height), QPointF(x, height - p.intensity2 * sI));
       }
    }    
  }
  else {
    // Draw raw peak data without colors and annotations, etc.
    vector<ms2::peak2> &ms2 = peak->ms2;
    painter->setPen(QPen(white));
    // Display only peaks within current range.
    for(size_t i = 0; i < ms2.size(); i++) {
      lcms::peak2 &p = ms2[i];
      if(p.mz2 < maxmz && p.mz2 > minmz) {
	double x = (p.mz2 - minmz) * smz;
	painter->drawLine(QPointF(x, height), QPointF(x, height - p.intensity2 * sI));
      }
    }
  }

}

// Add the MSMS View to the current scen, set width and height to fit
// into current graphics view.
MSMSWin::MSMSWin(lcms::analysis *a, QWidget *parent) : QDialog(parent) { 
  setupUi(this); 
  analysis = a;
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  graphicsView->setScene(&scene);
  int fw = graphicsView->frameWidth();
  view = new MSMSView(analysis, graphicsView->width() - 2*fw - 1, graphicsView->height() - 2*fw - 1, a->min_mz2, a->max_mz2);
  scene.addItem(view);
  initRangeSet = false;
}

// Set the axes to their defulat range and update edit boxes.
void MSMSWin::on_resetAxesButton_clicked(bool) {  
  view->setDefaultRange(); 
  minmzEdit->setText(QString::number(view->getMinmz()));
  maxIEdit->setText(QString::number(view->getMaxI()));
  maxmzEdit->setText(QString::number(view->getMaxmz()));
  view->update(); 
}

// Parse the text in each of the edit boxes and update the current
// range of the MS/MS view.
void MSMSWin::on_updateAxesButton_clicked(bool) {
  bool ok1 = false, ok2 = false, ok3 = false;
  double minmz = minmzEdit->text().toDouble(&ok1);
  double maxmz = maxmzEdit->text().toDouble(&ok2);
  double maxI = maxIEdit->text().toDouble(&ok3);
  if(ok1 && ok2 && ok3) {
    view->setRange(minmz, maxmz,maxI);
    view->update();
  }
}

// Adds peak to MS/MS display window.
void MSMSWin::setPeak(lcms::ms2peak *p, string descr) {
  if(p == NULL) return;

  // Add row to peak history list.
  QString displaydescr = descr.c_str();
  QString disp = displaydescr + ":(" + QString::number(p->retentionTime, 'f') + "," +  QString::number(p->mz, 'f') + ")";
    
  historyList->addItem(disp);

  peaks.push_back(p); // Update peak list.

  if(p->iso == ms2::isoHEAVY) view->setPeak(p, analysis->DB[ms2::dbHEAVY]);   // Display the new peak.
  else if(p->iso == ms2::isoMEDIUM) view->setPeak(p, analysis->DB[ms2::dbMEDIUM]);
  else view->setPeak(p, analysis->DB[ms2::dbLIGHT]);

  setChargeLabel(p->charge); // Change charge.
  setDBFragments(p); // Populate fragment list.
  peak_idx = peaks.size() - 1; // Update current peak index.

  // Update index of currently displayed fragment.
  if(p->results.size() == 0) frag_idx.push_back(-1); 
  else frag_idx.push_back(0);

  // Update selected item in history list.
  historyList->setCurrentRow(peak_idx);

  // Reinitalize range if necessary.
  if(initRangeSet == false) {
    if(p->ms2.size() > 0) {
      view->setDefaultRange();
      minmzEdit->setText(QString::number(view->getMinmz()));
      maxmzEdit->setText(QString::number(view->getMaxmz()));
      maxIEdit->setText(QString::number(view->getMaxI()));
      initRangeSet = true;
    }
  }

  view->update();
}


// Update peak if history list selection modified.
void MSMSWin::on_historyList_currentRowChanged(int row) { 
  if(row >= 0 && peaks.size() > 0 && row < (int)peaks.size()) {
    peak_idx = row;
    // Update current peak displayed and update charge.
    if(peaks[peak_idx]->iso == ms2::isoHEAVY) view->setPeak(peaks[peak_idx], analysis->DB[ms2::dbHEAVY]); 
    else if(peaks[peak_idx]->iso == ms2::isoMEDIUM) view->setPeak(peaks[peak_idx], analysis->DB[ms2::dbMEDIUM]);
    else view->setPeak(peaks[peak_idx], analysis->DB[ms2::dbLIGHT]); 

    // Set charge.
    setChargeLabel(peaks[peak_idx]->charge);

    // Change list of DB fragments.
    setDBFragments(peaks[peak_idx]);
    // Set current fragment using saved memory of fragments.
    if(frag_idx[peak_idx] >= 0) 
      fragmentList->setCurrentRow(frag_idx[peak_idx]);  
    else
      view->setFragment(NULL, NULL);
  }
}

// Clear history list. Reset MS/MS view state.
void MSMSWin::on_clearHistoryButton_clicked(bool) {  
  historyList->clear(); 
  fragmentList->clear();
  peaks.clear(); 
  frag_idx.clear();
  view->setPeak(NULL, NULL); 
  view->setFragment(NULL, NULL);
  initRangeSet = false;
  setChargeLabel(0);
}


void MSMSWin::on_fragmentList_currentRowChanged(int row) {
  if(row >= 0 && row < (int)peaks[peak_idx]->results.size() && frag_idx[peak_idx] >= 0) {
    view->setFragment(&peaks[peak_idx]->results[row].frag, &peaks[peak_idx]->results[row].mods);
    frag_idx[peak_idx] = row;
  }
}

// Update list of tryptic fragments.
void MSMSWin::setDBFragments(lcms::ms2peak *p) {
  fragmentList->clear();
  if(p == NULL) return;
  for(size_t i = 0; i < p->results.size(); i++) {
    QString score = QString::number(p->results[i].score);
    QString seqStr = analysis->get_seqS(&(p->results[i].frag), p->results[i].mods).c_str();
    QString meta = analysis->get_meta(&(p->results[i].frag)).substr(0,20).c_str();
    fragmentList->addItem(meta + ":" + seqStr + ":" + score);
  }
}
