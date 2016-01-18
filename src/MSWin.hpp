#ifndef __MSWIN_HPP__
#define __MSWIN_HPP__

#include <QtGui>
#include <QtWidgets>
#include <list>
#include <map>
#include <vector>

#include "ui_MSWin.h"

#include "dialogs.hpp"
#include "MSMSWin.hpp"
#include "analysis.hpp"

#include <iostream>
using namespace std;

// Code is organized to keep all GUI state in the view LCMSView.  The
// window MSWin only alters this state based on user input.

// Enums used as state information in LCMSView.
enum queryrect_t { NO_QRECT, FILTER_QRECT, NEIGHBOR_QRECT };

// Used to set selection mode.
enum selectmode_t { NO_SELECT, ZOOM, MS2_SELECT, XIC_SELECT, XIC_MS2_SELECT };

// retention time range 
struct time_range { double min_time, max_time; };

// Datatype for m/z range 
struct mz_range {  double min_mz, max_mz; };

// GraphicsView framework item for displaying LC-MS/MS data.
// This class contains most of the GUI's state information.
class LCMSView : public QObject, public QGraphicsItem {
  Q_OBJECT
public:
  // Constructor creates QImage canvas, initializes LCMS view state.
  LCMSView(lcms::analysis *a, int w, int h, QLabel *mzLabel_ = NULL, QLabel *rtLabel_ = NULL) {
    analysis = a; 
    width = w; height = h; 
    mzLabel = mzLabel_;  rtLabel = rtLabel_;

    condIndex = 0; // condition number, show 0th (first) condition
    repIndex = 0;  // replicate number, show 0th (first) replicate
    irunIndex = 0; // current irun to display, show 0th irun

    // Set range to original range specified in config.xml.
    d.min_mz = analysis->min_mz; d.max_mz = analysis->max_mz;
    t.min_time = analysis->min_time; t.max_time = analysis->max_time;

    terrain_colors(colors, 256);     // Initialize GUI color set.

    // Use a QImage canavs (faster to use setPixel).
    canvas = new QImage(width, height, QImage::Format_ARGB32_Premultiplied);

    // Initial GUI state parameters.
    selection = NULL;  
    select_mode = NO_SELECT;
    display_raw = false; display_graph = false;
    display_qrect = NO_QRECT;
    display_corr = false; display_obs = false;  display_XICs = false;
    display_isotopic = false;  
    display_group = false;
    display_IDs = true;

    peak_size = 1;
    ms2win = NULL; xicwin = NULL;
    setAcceptHoverEvents(true); // Allows display of mouse coordinates in corner.
    redraw();
  }

  ~LCMSView() { 
    delete canvas; 
    if(selection != NULL) delete selection;  
    if(ms2win != NULL) delete ms2win; 
    if(xicwin != NULL) delete xicwin;
  }

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *) {
    if(canvas != NULL) painter->drawImage(0,0, *canvas);
  }
  QRectF boundingRect() const { return QRectF(0, 0, width, height); }
  void setSize(int w, int h) {  
    width = w;  height = h;  
    delete canvas;
    canvas = new QImage(w,h, QImage::Format_ARGB32_Premultiplied);
    redraw();
  }
  // START Accessor functions for controlling the GUI state.
  QString getDescription() { 
    QString a = analysis->conditions[condIndex]->description.c_str(); 
    QString d = analysis->conditions[condIndex]->replicates[repIndex]->description.c_str(); 
    return a + ":" + d;
  }
  void setZoomMode() { select_mode = ZOOM; }
  void setMS2Select() { select_mode = MS2_SELECT; }
  void setXICSelect() { select_mode = XIC_SELECT; } 
  void setXICMS2Select() { select_mode = XIC_MS2_SELECT; }
  void zoomOut() { 
    if(zoom_mz.size() == 0) return;
    else {
      d = zoom_mz.front(); t = zoom_time.front(); 
      zoom_mz.pop_front(); 
      zoom_time.pop_front(); 
      redraw(); update();
    }
  }    
  void toggleRawData() { display_raw = !display_raw; redraw(); update(); }
  void toggleGraph() { display_graph = !display_graph; redraw(); update(); }
  void toggleGroup() { display_group = !display_group; redraw(); update(); }
  void incrementPeakSize() { peak_size++; redraw(); update(); }
  void decrementPeakSize() { if(peak_size > 1) { peak_size--; redraw(); update(); } }
  void cycleRect() { // Cycle through query rectangle display.
    switch(display_qrect) {
    case NO_QRECT: display_qrect = FILTER_QRECT; break;
    case FILTER_QRECT: display_qrect = NEIGHBOR_QRECT; break;
    case NEIGHBOR_QRECT: display_qrect = NO_QRECT; break;
    }
    redraw(); update();
  }

  void cycleCorrRect() { display_corr = !display_corr;  redraw(); update();  }
  void toggleXICs() { display_XICs = !display_XICs; redraw(); update(); }
  void toggleIDs() { display_IDs = !display_IDs; redraw(); update(); }
  void toggleIsotopic() { display_isotopic = !display_isotopic; redraw(); update(); }
  void toggleObs() { display_obs = !display_obs; redraw(); update(); }

  // Change or access the current condition, replicate set, or irun.
  void setCond(int newCondIndex) { condIndex = newCondIndex; }
  int getCond() { return condIndex; }
  void setRep(int newRepIndex) { repIndex = newRepIndex; }
  int getRep() { return repIndex; }
  void setIrun(int newIrunIndex) { irunIndex = newIrunIndex; }
  int getIrun() { return irunIndex; }

  // Change view position.
  void setmzrt(double mz1, double mz2, double rt1, double rt2) { 
    zoom_mz.push_front(d); zoom_time.push_front(t); // Save the zoom value.
    d.min_mz = mz1; d.max_mz = mz2; 
    t.min_time = rt1; t.max_time = rt2;
    redraw(); update();
  }
  // END accesor functions for controlling GUI state.

  void redraw(); 
protected:
  QLabel *mzLabel, *rtLabel; // pointers to label objects in MSWin updated as mouse moves.
  int condIndex, repIndex, irunIndex; // current condition, replicate, or irun to display

  lcms::analysis *analysis; // condition/replicate/irun LC-MS/MS data set

  mz_range d; time_range t; // Current m/z and time range.
  int width, height; // Current width and height of canvas/drawing area.
  list<mz_range> zoom_mz;  // Stack of ranges vales for different zoom settings for m/z dimension
  list<time_range> zoom_time; // and for the retention time dimension.
  vector<QColor> colors; // Coloring scheme for peak intensity.

  selectmode_t select_mode; // Switch to select mode.
  bool display_raw; // Toggles display of raw peak data.
  bool display_graph; // Toggles display of peak graph.
  bool display_XICs; // Flag turns on and off the display of XICs.
  bool display_isotopic; // Flag turns on and off display isotope/charge related boxes. 
  bool display_obs; // Flag for display of grouped XICs.
  bool display_group; // Display time range used to group across iruns.
  bool display_IDs; // Display search IDs.

  queryrect_t display_qrect; // Controls display of query rectangles.
  bool display_corr; // Controls display of alignment/grouping rectangles.
  
  int peak_size; // Size of drawn peak. 

  MSMSWin *ms2win; // non-modal dialog window for display MS/MS spectra.
  XICWin *xicwin; // non-modal dialog window for display XICs

  // Data drawing routines.
  QImage *canvas;
  void draw_XICs(lcms::irun *irun, QPainter *painter);
  void draw_isotope_groups(lcms::irun *irun, QPainter *painter);
  void draw_peaks(vector<lcms::peak *> &data, vector<lcms::ms2peak *> &data2, QPainter *painter);
  void draw_qrect(lcms::irun *irun, vector<lcms::peak *> &data, QPainter *painter);
  void draw_graph(lcms::irun *irun, QPainter *painter);
  void draw_grouped(QPainter *painter);
  void draw_grouped(QPainter *painter, vector<lcms::align_group *> &obs, QColor color);

  // Convert from image coordinates to a (m/z, retention time group).
  void from_canvas(double x, double y, double &rt, double &mz){      
    rt = x / double(width) * (t.max_time - t.min_time) + t.min_time;
    mz = y / double(height) * (d.max_mz - d.min_mz) + d.min_mz;
  }

  // Take an (m/z, time) group and convert to a canvas coordinate.
  void to_canvas(double retentionTime, double mz, int &x, int &y){
    double scaletime = double(width) / (t.max_time - t.min_time);
    double scalemz = double(height) / (d.max_mz - d.min_mz);
    x = (int)((retentionTime - t.min_time) * scaletime);
    y = (int)((mz - d.min_mz) * scalemz );      
  }

  // Selection rectangle, used for zooming into the data.
  QGraphicsRectItem *selection;
  // Mouse events used for zooming.
  void mousePressEvent( QGraphicsSceneMouseEvent *event) {
    if(selection == NULL && select_mode != NO_SELECT) {
      selection = new QGraphicsRectItem(event->pos().x(), event->pos().y(), 1,1, this);
      selection->setPen(QPen(Qt::white));
    }
  }
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event) {
    if(selection != NULL) {
      QRectF r = selection->rect();
      r.setRight(event->pos().x());
      r.setBottom(event->pos().y());
      selection->setRect(r);
    }
  }
  void mouseReleaseEvent( QGraphicsSceneMouseEvent *);

  void hoverMoveEvent(QGraphicsSceneHoverEvent *event) {
    // Capture mouse position, show corresponding (rt, m/z) position.
    double rt, mz;
    from_canvas(event->pos().x(), event->pos().y(), rt, mz);
    QString rtStr = QString::number(rt,'f',7);
    QString mzStr = QString::number(mz,'f',7);
    if(mzLabel != NULL) mzLabel->setText(mzStr);
    if(rtLabel != NULL) rtLabel->setText(rtStr);
  }
};


// Main window that manages the current set of replicate LC-MS/MS
// iruns.  This is just a wrapper around LCMSView which has most of
// the GUI state information.
class MSWin : public QMainWindow, private Ui::MSWin {
    Q_OBJECT
public:
  MSWin(QMainWindow *parent = 0);
  ~MSWin() { 
    if(view != NULL) delete view; 
    if(analysis != NULL) delete analysis;
  }

protected:
  QGraphicsScene scene;  // Scene that contains LCMSView  

  // Selected replicate, irun memory information.
  vector<int> selectedRepIndex;  // condIndex -> repIndex 
  vector< vector<int> > selectedIrunIndex; // (condIndex, repIndex) -> irunIndex 

  LCMSView *view; // Main LC-MS/MS view
  msconfig conf; // Global configuration information.
  lcms::analysis *analysis; // Current replicate set.

  // Mapping from tree widget (pointer) item to a specific isotope group.
  map<QTreeWidgetItem *, lcms::isotope_group*> item2group;
  // Mapping from atree widget to a specific XIC. 
  map<QTreeWidgetItem *, lcms::xic*> item2xic;
  
  void resizeEvent(QResizeEvent *) {
    if(view != NULL) {
      // Keep view the same size as window.
      int fw = graphicsView->frameWidth();
      view->setSize(graphicsView->width() - 2*fw-1, graphicsView->height() - 2*fw-1);
      graphicsView->setSceneRect(0,0,graphicsView->width() - 2*fw-1, graphicsView->height() - 2*fw-1);
    }
  }
  // Examines given directory. Collects subdirectories = conditions and mzXML files.
  // Keyboard capture.
  void keyReleaseEvent(QKeyEvent *);
  // Collect mzXML files in given directory.
  void collect_mzXML_files(QDir dir, vector<lcms::irun_file> &irunFiles);
  // Examines given directory for files with given suffix.
  void find_files(QString dir, vector<string> &base, vector<string> &files, QString suf);
  // Process iruns, divided among threads.
  void processIruns();
  // Multi-threaded irun/irun alignment.
  void alignIruns();
  // Loads fasta files in directory and builds fragment database, etc. 
  void loadDBs();
  // Utility function for displaying progress dialog.
  void progressDialog(QString message, QThread *thread) {
    vector<QThread *> threads; threads.push_back(thread);
    progressDialog(message, threads);
  }
  void progressDialog(QString message, vector<QThread *> threads);
  bool profileDialog(bool &ms1_profile, bool &ms2_profile);

  // Populates tree widgets.
  void populateIsotopeTree();
  void populateIsotopeGroups(QTreeWidgetItem *frags, vector<lcms::isop_seqgroup> &seq_groups);
  void populateXICTreeAlign();
  void populateXICTree();
  // Generates subtrees used by above functions.
  QTreeWidgetItem *accessionSubtree(vector<int> &idxs);
  QTreeWidgetItem *xicSubtree(string seq, int charge, vector<lcms::xic *> &xics);

  void closeEvent(QCloseEvent *event) {
    if(view != NULL) { delete view; view = NULL; } 
    if(analysis != NULL) { delete analysis; analysis = NULL; }
    QMainWindow::closeEvent(event);
  }
    

private slots:
  void on_action_Open_triggered();
  void on_action_Recalibrate_triggered();
  void on_action_Quit_triggered() {
    if(view != NULL) delete view; if(analysis != NULL) delete analysis;
    exit(0);
  }
  void on_actionSave_TXT_triggered();
  void on_actionSave_SILAC_TXT_triggered();
  void on_actionSave_Align_XICs_TXT_triggered();

  void on_isotopeTree_itemDoubleClicked(QTreeWidgetItem * item, int column);
  void on_xicTree_itemDoubleClicked(QTreeWidgetItem * item, int column);

  void on_listCond_currentRowChanged(int row);
  void on_listRep_currentRowChanged(int row);
  void on_listIrun_currentRowChanged(int row);
};

#endif // __MSWIN_HPP__
