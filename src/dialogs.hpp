#ifndef __DIALOGS_HPP__
#define __DIALOGS_HPP__

#include <QtGui>
#include <QtWidgets>
#include <list>
#include <vector>

#include "ui_MSConfigWin.h"
#include "ui_XICWin.h"
#include "ui_MassWin.h"
#include "ui_AAWin.h"
#include "ui_RecalWin.h"

#include "analysis.hpp"

#include <iostream>
using namespace std;


// Generate colors that have blue at the lowest level and go through
// green up to yellow and brown to display log(intensity) values.  
inline void terrain_colors(vector<QColor> &colors, int n){
  QColor C;
  double h = 258, s = 1, v = 1;   // Generate the blue range.
  for(int i = 0; i < n - n/3 - n/3; i++){
    h -= (258 - 186) / (double(n - n/3 - n/3));
    C.setHsv((int)h, int(255*s), int(255*v));
    colors.push_back(C);
  }
  h = 138; s = 1; v = 1;   // Generate the green range.
  for(int j = 0; j < n/3; j++){
    h -= (138 - 66) / (double(n/3));
    C.setHsv((int)h, int(255*s), int(255*v));
    colors.push_back(C);
  }
  h = 60; s = 1; v = 1;   // Generate the yellow/brown range.
  for(int k = 0; k < n/3; k++){
    h -= (60 - 36) / (double(n/3)); s -= (1 - 0.3) / (double(n/3));
    C.setHsv((int)h, int(255*s), int(255*v));
    colors.push_back(C);
  }
}

// Provides a simple calculator for computing mass shift.
// Note: Mass spectrometrists cannot be trusted to enter floating point values!!!
class MassWin : public QDialog, private Ui::MassWin {
  Q_OBJECT
public:
  MassWin(QWidget *parent = 0) : QDialog(parent) { 
    setupUi(this); _13C = _15N = _2H = _18O =  _17O = C = H = O = N = S = P = 0;
  }
  double getMassDouble() {
    double calcMass = C * mass::C + H * mass::H + O * mass::O + S * mass::S + P * mass::P + N * mass::N;
    using namespace mass;
    double isoMass = _13C * (delta13C + mass::C) + _15N * (delta15N + mass::N) + 
      _2H * (delta2H + mass::H) + _18O * (delta18O + mass::O) + _17O * (delta17O + mass::O);
    return isoMass + calcMass;
  }
  QString getMassStr() { return QString::number(getMassDouble(), 'f', 10); }
protected:
  double C, H, O, N, S, P;
  double _13C, _15N, _2H, _18O, _17O; // isotopes
  void updateMass() {labelMass->setText(getMassStr()); }
private slots:
  void on_spin13C_valueChanged(int i) { _13C = i; updateMass(); } void on_spin15N_valueChanged(int i) { _15N = i; updateMass(); }
  void on_spin2H_valueChanged (int i) { _2H  = i; updateMass(); } void on_spin18O_valueChanged(int i) { _18O = i; updateMass(); }
  void on_spin17O_valueChanged(int i) { _17O = i; updateMass(); } void on_spinC_valueChanged(int i) { C = i; updateMass(); }
  void on_spinH_valueChanged(int i) { H = i; updateMass(); } void on_spinO_valueChanged(int i) { O = i; updateMass(); }
  void on_spinN_valueChanged(int i) { N = i; updateMass(); } void on_spinS_valueChanged(int i) { S = i; updateMass(); }
  void on_spinP_valueChanged(int i) { P = i; updateMass(); }
};

// Again, mass spectrometrists should not be trusted with the IUPAC amino acid codes as well.
class AAWin : public QDialog, private Ui::AAWin {
  Q_OBJECT
public:
  AAWin(QWidget *parent = 0) : QDialog(parent) { 
    setupUi(this);
    listAAs->addItem("A = alanine");       listAAs->addItem("C = cysteine");
    listAAs->addItem("D = aspartic acid"); listAAs->addItem("E = glutamic acid");
    listAAs->addItem("F = phenylalanine"); listAAs->addItem("G = glycine");
    listAAs->addItem("H = histidine");     listAAs->addItem("I = isoleucine");
    listAAs->addItem("K = lysine");        listAAs->addItem("L = leucine");
    listAAs->addItem("M = methionine");    listAAs->addItem("N = asparagine");
    listAAs->addItem("P = proline");       listAAs->addItem("Q = glutamine");
    listAAs->addItem("R =  arginine");     listAAs->addItem("S = serine");
    listAAs->addItem("T = threonine");     listAAs->addItem("V = valine");
    listAAs->addItem("W = tryptophan");    listAAs->addItem("Y = tyrosine");
    listAAs->addItem("[ = cleavage N-terminus");
    listAAs->addItem("] = cleavage C-terminus");
    listAAs->addItem("# = protein N-terminus");
    listAAs->addItem("* = protein C-terminus");    
  }
  QString getSelAAs() {
    QString sel;
    if(listAAs->item(0)->isSelected())  sel += 'A'; if(listAAs->item(1)->isSelected())  sel += 'C';
    if(listAAs->item(2)->isSelected())  sel += 'D'; if(listAAs->item(3)->isSelected())  sel += 'E';
    if(listAAs->item(4)->isSelected())  sel += 'F'; if(listAAs->item(5)->isSelected())  sel += 'G';
    if(listAAs->item(6)->isSelected())  sel += 'H'; if(listAAs->item(7)->isSelected())  sel += 'I';
    if(listAAs->item(8)->isSelected())  sel += 'K'; if(listAAs->item(9)->isSelected())  sel += 'L';
    if(listAAs->item(10)->isSelected()) sel += 'M'; if(listAAs->item(11)->isSelected()) sel += 'N';
    if(listAAs->item(12)->isSelected()) sel += 'P'; if(listAAs->item(13)->isSelected()) sel += 'Q';
    if(listAAs->item(14)->isSelected()) sel += 'R'; if(listAAs->item(15)->isSelected()) sel += 'S';
    if(listAAs->item(16)->isSelected()) sel += 'T'; if(listAAs->item(17)->isSelected()) sel += 'V';
    if(listAAs->item(18)->isSelected()) sel += 'W'; if(listAAs->item(19)->isSelected()) sel += 'Y';
    if(listAAs->item(20)->isSelected()) sel += '['; if(listAAs->item(21)->isSelected()) sel += ']';
    if(listAAs->item(22)->isSelected()) sel += '#'; if(listAAs->item(23)->isSelected()) sel += '*';
    return sel;
  }
};


// Configuration dialog displayed prior to loading mzXML data files.
class MSConfigWin : public QDialog, private Ui::MSConfigWin {
  Q_OBJECT
public:
  MSConfigWin(msconfig &conf_, QWidget *parent = 0);
protected:
  msconfig &conf;
  bool filter_param_changed_saved;
  void PopulateFixedMods();
  void PopulateVarMods();
  void PopulateIsotopes();
  QString IsotopeText(isotope_t &iso);
  QString FixedModText(fxmod_t &fx);
  QString VarModText(varmod_t &vr);
private slots:
  // Button handlers below. 
  void on_cancelButton_clicked(bool) {
    conf.filter_param_changed = filter_param_changed_saved;
    reject(); 
  }
  void on_loadButton_clicked(bool);

  void on_addIsotopeButton_clicked(bool);  void on_removeIsotopeButton_clicked(bool);
  void on_addFixedButton_clicked(bool);    void on_removeFixedButton_clicked(bool);
  void on_addFragShiftButton_clicked(bool);void on_removeFragShiftButton_clicked(bool);
  void on_addVarButton_clicked(bool);      void on_removeVarButton_clicked(bool);

  // Mass calculation and amino acid selection buttons.
  void on_calcHeavyButton_clicked(bool) {
    MassWin calc(this); if(calc.exec() == QDialog::Rejected) return; 
    editIsoShiftHeavy->setText(calc.getMassStr());  
  }
  void on_calcMediumButton_clicked(bool) {
    MassWin calc(this); if(calc.exec() == QDialog::Rejected) return; 
    editIsoShiftMedium->setText(calc.getMassStr());
  }
  void on_calcFixedButton_clicked(bool) {
    MassWin calc(this); if(calc.exec() == QDialog::Rejected) return; 
    editFixedShift->setText(calc.getMassStr());
  }
  void on_calcPMassButton_clicked(bool) {
    MassWin calc(this); if(calc.exec() == QDialog::Rejected) return; 
    editDeltaPMass->setText(calc.getMassStr());
  }
  void on_isoAAsButton_clicked(bool) {
    AAWin aawin(this); if(aawin.exec() == QDialog::Rejected) return; 
    editIsotopeAAs->setText(aawin.getSelAAs());
  }
  void on_calcFragMassButton_clicked(bool) {
    MassWin calc(this); if(calc.exec() == QDialog::Rejected) return; 
    editFragShift->setText(calc.getMassStr());
  }
  // Amino acid list entry buttons
  void on_fixedAAsButton_clicked(bool) {
    AAWin aawin(this); if(aawin.exec() == QDialog::Rejected) return; 
    editFixedAAs->setText(aawin.getSelAAs());
  }
  void on_fragAAButton_clicked(bool){
    AAWin aawin(this); if(aawin.exec() == QDialog::Rejected) return; 
    editFragAAs->setText(aawin.getSelAAs());
  }
};


class RecalWin : public QDialog, private Ui::RecalWin {
  Q_OBJECT
public:
  RecalWin(QWidget *parent = 0) : QDialog(parent) { recalMzFlag = recalTimeFlag = false; setupUi(this); }
  double getMzWin() { return spinMzWindow->value(); } 
  double getTimeWin() { return spinTimeWindow->value(); }
  bool recalMzClicked() { return recalMzFlag; }  bool recalTimeClicked() { return recalTimeFlag; }
private:
  bool recalMzFlag, recalTimeFlag;
private slots:
  // Button handlers below. 
  void on_mzRecalButton_clicked(bool) { recalMzFlag = true; accept(); }
  void on_timeRecalButton_clicked(bool) { recalTimeFlag = true; accept(); }
  void on_cancelButton_clicked(bool) { reject(); }
};


// Draws an MS/MS fragmentation spectrum. 
class XICView : public QObject, public QGraphicsItem {
  Q_OBJECT
public:
  XICView(msconfig &conf_, int w, int h) : conf(conf_) {  
    width = w; height = h; xic = NULL; 
    setAcceptHoverEvents(true); // Allows mouse position display in (rt, I) coordinates.
    positionText = NULL;
    select_mode = false;
    selection = NULL;
    terrain_colors(colors, 256);
  }
  ~XICView() { if(positionText != NULL) delete positionText; }
  void paint(QPainter *, const QStyleOptionGraphicsItem *, QWidget *);
  QRectF boundingRect() const { return QRectF(0, 0, width, height); }  
  void setSize(int w, int h) {  width = w;  height = h;  }
  // Sets the m/z and intensity ranges to maximum and minimum values.
  void setDefaultRange() {
    if(xic == NULL) return;

    maxrt_stack.clear(); minrt_stack.clear(); maxI_stack.clear();

    maxI = -(numeric_limits<double>::max());
    minrt = numeric_limits<double>::max(); maxrt = -(numeric_limits<double>::max());
    for(size_t i = 0; i < xic->chrom.size(); i++) {
      double rt = xic->chrom[i].retentionTime;
      double I = xic->chrom[i].intensity;
      if(rt < minrt) minrt = rt; if(rt > maxrt) maxrt = rt;
      if(I > maxI) maxI = I;
    }
  }
  // Modify and access the range data.
  void setRange(double minrt_, double maxrt_, double maxI_) { minrt = minrt_; maxrt = maxrt_; maxI = maxI_;  }
  double getMinrt() { return minrt; }
  double getMaxrt() { return maxrt; }
  double getMaxI() { return maxI; }
  // Update the current XIC.
  void setXIC(lcms::xic *x) { 
    if(x == NULL) return;
    xic = x; 
    update(); 
  }
  void setZoomMode() { select_mode = true;  }
  void zoomOut() {
    if(maxrt_stack.size() >= 1){
      maxrt = maxrt_stack.front(); maxrt_stack.pop_front();
      minrt = minrt_stack.front(); minrt_stack.pop_front();
      maxI = maxI_stack.front(); maxI_stack.pop_front();
      update();
    }
  }
protected:
  msconfig &conf;
  QGraphicsTextItem *positionText;
  int width, height;   // view width & height
  double maxI;   // y-axis intensity range
  double minrt, maxrt; // x-axis m/z range
  lcms::xic *xic; 
  list<double> maxrt_stack;
  list<double> minrt_stack;
  list<double> maxI_stack;
  
  QGraphicsRectItem *selection;
  bool select_mode;

  vector<QColor> colors; // Coloring scheme for peak intensity.
  
  void hoverMoveEvent(QGraphicsSceneHoverEvent *event){ 
    // Translate mouse coordinates in MS/MS view into (m/z, intensity).
    double I = double(height - event->pos().y()) / double(height) * maxI;
    double rt = event->pos().x() / double(width) * (maxrt - minrt) + minrt;
    // Update position string on display.
    QString posStr = QString::number(rt,'f',6) + ", " + QString::number(I,'g',6);
    if(positionText == NULL) { 
      positionText = new QGraphicsTextItem(posStr, this);
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
      maxrt_stack.push_front(maxrt);
      minrt_stack.push_front(minrt);
      maxI_stack.push_front(maxI);

      // Convert into (rt, I) coordinates. 
      double rtmin_new = rect.left() / double(width) * (maxrt - minrt) + minrt;
      double maxrt_new = rect.right() / double(width) * (maxrt - minrt) + minrt;
      double maxI_new = double(height - rect.top()) / double(height) * maxI;

      minrt = rtmin_new;
      maxrt = maxrt_new;
      maxI = maxI_new;

      // Update the current zoom region.
      update();
      
      delete selection;
      selection = NULL; select_mode = false;
    }
  }
};

// Main dialog manages current XIC view.
class XICWin : public QDialog, private Ui::XICWin {
  Q_OBJECT
public:
  XICWin(msconfig &conf_, QWidget *parent = 0);
  ~XICWin() { if(view != NULL) delete view; }
  void setXIC(lcms::xic *x, string descr) { 
    QString displaydescr = descr.c_str();
    QString disp = displaydescr + ":" + QString::number(x->quant, 'f', 0) +
      ":(" + QString::number(x->retentionTime, 'f', 1) + "," + QString::number(x->mz, 'f') + "):" + 
      "+" + QString::number(x->charge);
    historyList->addItem(disp);
    xics.push_back(x);
    dataset.push_back(descr);
    view->setXIC(x);
    if(initRangeSet == false) {
      view->setDefaultRange();
      minrtEdit->setText(QString::number(view->getMinrt()));
      maxrtEdit->setText(QString::number(view->getMaxrt()));
      maxIEdit->setText(QString::number(view->getMaxI()));
      initRangeSet = true;
    }
    labelCharge->setText("+" + QString::number(x->charge));
    // update labels, etc. 
    view->update();
  }

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
      }
    }
  }

protected:
  msconfig &conf;
  XICView *view; // XIC
  // NOTE: In the list item in window, display precusor: (m/z, retention time, intensity). 
  vector<lcms::xic *> xics; // List of previously displayed MS/MS peaks.
  vector<string> dataset; // Descriptions of data sets.
  bool initRangeSet; // Set to true if the initial MS/MS scan range has been set.
  QGraphicsScene scene;
private slots:
  void on_historyList_currentRowChanged(int);
  void on_clearHistoryButton_clicked(bool);
  void on_resetAxesButton_clicked(bool);
  void on_updateAxesButton_clicked(bool);
};


#endif
