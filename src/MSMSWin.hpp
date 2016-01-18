#ifndef __MSMSWin_HPP__
#define __MSMSWin_HPP__

#include <QtGui>
#include <QtWidgets>
#include <list>
#include <vector>

#include "ui_MSMSWin.h"
#include "analysis.hpp"

// Draws an MS/MS fragmentation spectrum. 
class MSMSView : public QObject, public QGraphicsItem {
  Q_OBJECT
public:
  MSMSView(lcms::analysis *a_ ,int w, int h, float m_mz2, float x_mz2) { 
    analysis = a_;
    min_mz2 = m_mz2; max_mz2 = x_mz2;
    ms2db = NULL;
    width = w; height = h; 
    peak = NULL; 
    positionText = NULL;
    selection = NULL;
    select_mode = false;
    setAcceptHoverEvents(true); // Allows mouse position display in (mz, I) coordinates.
  }
  ~MSMSView() { if(positionText != NULL) delete positionText; }

  void paint(QPainter *, const QStyleOptionGraphicsItem *, QWidget *);
  QRectF boundingRect() const { return QRectF(0, 0, width, height); }  
  void setSize(int w, int h) {  width = w;  height = h;  }

  // Sets the m/z and intensity ranges to maximum and minimum values.
  void setDefaultRange() {
    if(peak == NULL) return;
    maxI = -(numeric_limits<double>::max());
    minmz = numeric_limits<double>::max(); maxmz = -(numeric_limits<double>::max());
    vector<ms2::peak2> &ms2 = peak->ms2;
    for(size_t i = 0; i < ms2.size(); i++) {
      double mz = ms2[i].mz2;
      double I = ms2[i].intensity2;
      if(mz < minmz) minmz = mz;
      if(mz > maxmz) maxmz = mz;
      if(I > maxI) maxI = I;
    }
    maxI = maxI * 1.25;
  }

  // Modify and access the range data.
  void setRange(double minmz_, double maxmz_, double maxI_) { minmz = minmz_; maxmz = maxmz_; maxI = maxI_;  }
  double getMinmz() { return minmz; }
  double getMaxmz() { return maxmz; }
  double getMaxI() { return maxI; }

  void setFragment(ms2::fragment *f, ms2::mods_t *mod_) { fragment = f; mod = mod_;  update(); }

  // Update the current peak.
  void setPeak(lcms::ms2peak *p, ms2::fragmentDB *db) { 
    peak = p; 
    ms2db = db;
    if(peak == NULL) { maxmz_stack.clear(); minmz_stack.clear();  maxI_stack.clear(); }
    update(); 
  }
  void setZoomMode() { select_mode = true; }
  void zoomOut() {
    if(maxmz_stack.size() >= 1){
      maxmz = maxmz_stack.front(); maxmz_stack.pop_front();
      minmz = minmz_stack.front(); minmz_stack.pop_front();
      maxI = maxI_stack.front(); maxI_stack.pop_front();
      update();
    }
  }

protected:
  lcms::analysis *analysis;
  QGraphicsTextItem *positionText;
  int width, height;   // view width & height
  double maxI;   // y-axis intensity range
  double minmz, maxmz; // x-axis m/z range

  lcms::ms2peak *peak;  // Current MS/MS peak.
  ms2::fragmentDB *ms2db; // Global fragment database.
  ms2::fragment *fragment; // Current MS/MS theoretical fragment to display
  ms2::mods_t *mod; // Any modifications applied to the current fragment.
  float min_mz2, max_mz2;

  list<double> maxmz_stack, minmz_stack, maxI_stack; // Saved zoom information.
  
  QGraphicsRectItem *selection;
  bool select_mode;

  void hoverMoveEvent(QGraphicsSceneHoverEvent *event){ 
    // Translate mouse coordinates in MS/MS view into (m/z, intensity).
    double I = 0, mz = 0;
    if(peak != NULL) {
      // Compute value based on mouse position.
      if(peak->ms2.size() > 0) {
	I = (height - event->pos().y()) / double(height) * maxI;
	mz = event->pos().x() / double(width) * (maxmz - minmz) + minmz;
      }      
    }
    // Update position string on display.
    QString posStr = QString::number(mz,'f',6) + ", " + QString::number(I,'g',6);
    
    if(positionText == NULL) { 
      // Allocate position text object, if not allocated.
      positionText = new QGraphicsTextItem(posStr, this);
      QColor white(255,255,255);
      positionText->setDefaultTextColor(white);
      positionText->setPos(0,0);
    }
    positionText->setPlainText(posStr);
  }

  void mousePressEvent( QGraphicsSceneMouseEvent *event) {
    if(selection == NULL && select_mode) {
      selection = new QGraphicsRectItem(event->pos().x(), event->pos().y(), 1, height - event->pos().y(), this);
      selection->setPen(QPen(Qt::red));
    }
  }

  void mouseMoveEvent( QGraphicsSceneMouseEvent *event) {
    if(selection != NULL) {
      QRectF r = selection->rect();
      r.setRight(event->pos().x());
      selection->setRect(r);
    }
  }
  void mouseReleaseEvent( QGraphicsSceneMouseEvent *) {
    if(selection != NULL) {
      scene()->removeItem(selection);
      QRectF rect = selection->rect().normalized();

      // Save the zoom value.
      maxmz_stack.push_front(maxmz);
      minmz_stack.push_front(minmz);
      maxI_stack.push_front(maxI);

      // Convert into (mz, I) coordinates. 
      double mzmin_new = rect.left() / double(width) * (maxmz - minmz) + minmz;
      double maxmz_new = rect.right() / double(width) * (maxmz - minmz) + minmz;
      double maxI_new = double(height - rect.top()) / double(height) * maxI;

      minmz = mzmin_new;
      maxmz = maxmz_new;
      maxI = maxI_new;

      // Update the current zoom region.
      update();
      
      delete selection;
      selection = NULL; select_mode = false;
    }
  }

};


// Main dialog that manages current MS/MS scan.
class MSMSWin : public QDialog, private Ui::MSMSWin {
  Q_OBJECT
public:
  MSMSWin(lcms::analysis *a, QWidget *parent = 0);
  ~MSMSWin() { if(view != NULL) delete view; }

  void setPeak(lcms::ms2peak *p, string descr);

  void resizeEvent(QResizeEvent *) {
    int fw = graphicsView->frameWidth();
    view->setSize(graphicsView->width() - 2*fw-1, graphicsView->height() - 2*fw-1);
    graphicsView->setSceneRect(0,0,graphicsView->width() - 2*fw-1, graphicsView->height() - 2*fw-1);
  }

  void keyReleaseEvent(QKeyEvent *event) {
    if(view != NULL) {
      switch(event->key()) {
      case Qt::Key_Z: view->setZoomMode(); break;
      case Qt::Key_X: view->zoomOut(); break;
      default: break;
      }
    }
  }

protected:
  // Update label text that displays charge information.
  void setChargeLabel(int charge) {
    if(charge >= 0) chargeLabel->setText("Precursor charge: +" + QString::number(charge));
    else chargeLabel->setText("Precursor charge: " + QString::number(charge));
  }

  MSMSView *view; // MS/MS scan.
  lcms::analysis *analysis;
  vector<lcms::ms2peak *> peaks; // List of previously displayed MS/MS peaks.
  // NOTE: peaks is in one-to-one correspondence with historyList rows.
  int peak_idx; 

  vector<int> frag_idx;

  bool initRangeSet; // Set to true if the initial MS/MS scan range has been set.
  QGraphicsScene scene;

  // Update list of assign DB fragments.  
  void setDBFragments(lcms::ms2peak *p);

private slots:
  void on_resetAxesButton_clicked(bool);
  void on_updateAxesButton_clicked(bool);
  void on_historyList_currentRowChanged(int);
  void on_fragmentList_currentRowChanged(int);
  void on_clearHistoryButton_clicked(bool);
};



#endif
