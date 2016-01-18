#ifndef __THERMO_HPP__
#define __THERMO_HPP__

#include <QtGui>

#include "ui_thermo2pview.h"

class ThermoWin : public QMainWindow, private Ui::ThermoWin {
  Q_OBJECT
public:
  ThermoWin(QMainWindow *parent = 0) ;
  ~ThermoWin() {  
  }
protected:
private slots:
   void on_testButton_clicked(bool);
};

#endif
