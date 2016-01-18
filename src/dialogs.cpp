#include <iostream>
#include <iomanip>
#include <fstream>

#include <algorithm>
#include <map>
#include <ctype.h>

#include "dialogs.hpp"

#include <Qt>

// TODO: Add more checks on the validity of the parameters entered!!!

XICWin::XICWin(msconfig &conf_, QWidget *parent) : QDialog(parent), conf(conf_) { 
  setupUi(this); 
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  graphicsView->setScene(&scene);
  // Adjust for fram width. Create view.
  int fw = graphicsView->frameWidth();
  view = new XICView(conf, graphicsView->width() - 2*fw - 1, graphicsView->height() - 2*fw - 1);
  scene.addItem(view);
  initRangeSet = false;
}


void XICView::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *) {
  if(xic == NULL) return;

  QColor black(0,0,0);
  // Use black background for contrast.
  painter->setPen(QPen(black));
  painter->setBrush(QBrush(black));
  painter->drawRect(0,0, width, height);

  QColor green(0,255,0); painter->setPen(QPen(green));

  // Display XIC.
  double srt = width / (maxrt - minrt), sI = height / maxI;
  double logImin = conf.logImin_disp, logImax = conf.logImax_disp;
  
  // Display only points within current range.
  for(size_t i = 0; i < xic->chrom.size(); i++) {
    lcms::peak &p = xic->chrom[i];
    if(p.retentionTime < maxrt && p.retentionTime > minrt) {
      double x = (p.retentionTime - minrt) * srt;

      // Chose color to display.
      double logI = log10(p.intensity);
      int level = 0;
      if(logI > logImax) level = 255;
      else if( logI < logImin ) level = 0;
      else level = (int)((logI - logImin) / (logImax - logImin) * 256.0);

      QPen peak_color(colors[level]);
      painter->setPen(peak_color);
      painter->drawLine(QPointF(x, height), QPointF(x, height - p.intensity * sI));
    }
  }
  
}

// Update drawn XIC if history list selection modified.
void XICWin::on_historyList_currentRowChanged(int row) { 
  if(row >= 0 && xics.size() > 0 && row < (int)xics.size()) {
    view->setXIC(xics[row]); 
    labelCharge->setText("+" + QString::number(xics[row]->charge));
  }
}

// Clear XIC history list. 
void XICWin::on_clearHistoryButton_clicked(bool) {  
  historyList->clear(); 
  xics.clear(); 
  dataset.clear();
  labelCharge->setText("+0");
  view->setXIC(NULL); 
  initRangeSet = false;
}

// Set the axes to their defulat range and update edit boxes.
void XICWin::on_resetAxesButton_clicked(bool) {  
  view->setDefaultRange(); 
  minrtEdit->setText(QString::number(view->getMinrt()));
  maxrtEdit->setText(QString::number(view->getMaxrt()));
  maxIEdit->setText(QString::number(view->getMaxI()));
  view->update(); 
}

// Parse the text in each of the edit boxes and update the current
// range of the MS/MS view.
void XICWin::on_updateAxesButton_clicked(bool) {
  bool ok1 = false, ok2 = false, ok3 = false;
  double minrt = minrtEdit->text().toDouble(&ok1);
  double maxrt = maxrtEdit->text().toDouble(&ok2);
  double maxI = maxIEdit->text().toDouble(&ok3);
  if(ok1 && ok2 && ok3) {
    view->setRange(minrt, maxrt,maxI);
    view->update();
  }
}


// ----------------------------------------------------------
// ----------------------------------------------------------
// ----------------------------------------------------------

MSConfigWin::MSConfigWin(msconfig &conf_, QWidget *parent) : QDialog(parent), conf(conf_) {
  setupUi(this);

  filter_param_changed_saved = conf.filter_param_changed;
  conf.filter_param_changed = false;

  comboQuantMode->setCurrentIndex((int)conf.quantmode);
  spinThreads->setValue(conf.threads);

  // Set filtering parameters.
  editFilterdmz->setText(QString::number(conf.dmz));
  //editFilterdtime->setText(QString::number(conf.dtime));
  //editPeakThresh->setText(QString::number(conf.peakthresh));
  editHistThresh->setText(QString::number(conf.hist_log10I_thresh));

  // Set XIC parameters.
  //editXICdtime->setText(QString::number(conf.peak_dtime));
  editXICdmz->setText(QString::number(conf.peak_dmz));
  editXICmax->setText(QString::number(conf.max_xiclen));
  editXICWidth->setText(QString::number(conf.xic_width));
  editIsoTol->setText(QString::number(conf.isotolerence));
  spinXICPeakMissed->setValue(conf.xic_peaks_missed);
  spinMinXICPeaks->setValue(conf.min_xic_peaks);

  // Populate grouping parameters.
  if(conf.skip_normalize_label_free) checkSkipNormalizeLabelFree->setCheckState(Qt::Checked);
  else checkSkipNormalizeLabelFree->setCheckState(Qt::Unchecked);

  editAligndtime->setText(QString::number(conf.align_dtime));
  editMinConds->setText(QString::number(conf.min_conds));
  editMinReps->setText(QString::number(conf.min_reps));
  editMinIruns->setText(QString::number(conf.min_iruns));
  editGroupdtime->setText(QString::number(conf.group_dtime));

  if(conf.align_translation) checkAlignTranslation->setCheckState(Qt::Checked);
  else checkAlignTranslation->setCheckState(Qt::Unchecked);

  if(conf.align_nonlinear) checkAlignNonlinear->setCheckState(Qt::Checked);
  else checkAlignNonlinear->setCheckState(Qt::Unchecked);

  // Populate flags used to control memory usage in data processing.
  if(conf.run_algorithms) checkRunAlgo->setCheckState(Qt::Checked);
  else checkRunAlgo->setCheckState(Qt::Unchecked);

  if(conf.keep_raw) checkRawPeaks->setCheckState(Qt::Checked);
  else checkRawPeaks->setCheckState(Qt::Unchecked) ;

  if(conf.keep_filtered) checkFilteredPeaks->setCheckState(Qt::Checked);
  else checkFilteredPeaks->setCheckState(Qt::Unchecked);

  if(conf.load_ms2) checkLoadMS2->setCheckState(Qt::Checked);  
  else checkLoadMS2->setCheckState(Qt::Unchecked);  

  if(conf.keep_ms2) checkKeepMS2->setCheckState(Qt::Checked);
  else checkKeepMS2->setCheckState(Qt::Unchecked);

  // Populate MS/MS database search parameters.
  editMSTol->setText(QString::number(conf.mstol));
  editMS2Tol->setText(QString::number(conf.ms2tol));
  if(conf.ms2tol_usedaltons) comboMS2Units->setCurrentIndex(1); // daltons selected
  else comboMS2Units->setCurrentIndex(0); // ppms selected 

  editMissedCleaves->setText(QString::number(conf.missed_cleaves));

  // see indices in enzyme_t (msconfig.hpp)
  comboEnzyme->addItem("Trypsin");  // 0
  comboEnzyme->addItem("Lys-C"); // 1
  comboEnzyme->addItem("Glu-C"); // 2
  comboEnzyme->addItem("Arg-C"); // 3
  comboEnzyme->addItem("Trypsin/P"); // 4
  comboEnzyme->addItem("Lys-N"); // 5
  comboEnzyme->addItem("Chymotrypsin"); // 6
  comboEnzyme->addItem("Asp-N"); // 7

  comboEnzyme->setCurrentIndex(conf.enzyme_id);

  if(conf.fixed_15N) checkFixed15N->setCheckState(Qt::Checked);
  else checkFixed15N->setCheckState(Qt::Unchecked);

  editFDR->setText(QString::number(conf.ms2_fdr));

  PopulateFixedMods(); PopulateVarMods(); PopulateIsotopes();

  // Populate MS/MS mod related parameters
  editMaxModAA->setText(QString::number(conf.max_mod_aa));

  // Populate isotope related parameters.
  spinMaxLabel->setValue(conf.max_label_count);

  // Two metabolic labeling types.
  if(conf.isotope_15N) checkIsotope15N->setCheckState(Qt::Checked); 
  else checkIsotope15N->setCheckState(Qt::Unchecked);

  if(conf.isotope_by_search) checkUngroupedIsotopes->setCheckState(Qt::Checked);
  else checkUngroupedIsotopes->setCheckState(Qt::Unchecked);

  if(conf.skip_normalize_ratios) checkSkipNormalizeIsotope->setCheckState(Qt::Checked);
  else checkSkipNormalizeIsotope->setCheckState(Qt::Unchecked);

  if(conf.acetyl_stoich_plus5) checkAcetylPlus5->setCheckState(Qt::Checked);
  else checkAcetylPlus5->setCheckState(Qt::Unchecked);
}

void MSConfigWin::on_loadButton_clicked(bool) {
  conf.quantmode = (quantmode_t)comboQuantMode->currentIndex();

  // Number of CPU threads to use.
  conf.threads = spinThreads->value();

  // Collect filtering parameters.
  double dmz_new = editFilterdmz->text().toDouble();
  //double dtime_new = editFilterdtime->text().toDouble();
  //int peakthresh_new = editPeakThresh->text().toInt();
  double hist_log10I_thresh_new = editHistThresh->text().toDouble();

  // Keep track of whether any of these paramters changed.
  if(conf.dmz != dmz_new) conf.filter_param_changed = true; 
  //if(conf.dtime != dtime_new) conf.filter_param_changed = true;
  //if(conf.peakthresh != peakthresh_new) conf.filter_param_changed = true; 
  if(conf.hist_log10I_thresh != hist_log10I_thresh_new) conf.filter_param_changed = true; 

  conf.dmz = dmz_new;
  //conf.dtime = dtime_new; 
  //conf.peakthresh = peakthresh_new;
  conf.hist_log10I_thresh = hist_log10I_thresh_new;

  // Validate filtering parameters.
  //if( conf.dtime < 0 || conf.dmz < 0 || conf.peakthresh < 0) {  return;  }
  if( conf.dmz < 0 ) {  return;  }
  
  // Collect XIC parameters.
  //conf.peak_dtime = editXICdtime->text().toDouble();
  conf.peak_dmz = editXICdmz->text().toDouble();
  conf.xic_width = editXICWidth->text().toDouble();
  conf.max_xiclen = editXICmax->text().toDouble();
  conf.isotolerence = editIsoTol->text().toDouble();
  conf.xic_peaks_missed = spinXICPeakMissed->value();
  conf.min_xic_peaks = spinMinXICPeaks->value();

  conf.align_translation = checkAlignTranslation->checkState() == Qt::Checked;
  conf.align_nonlinear = checkAlignNonlinear->checkState() == Qt::Checked;

  // Validate XIC parameters.
  if(/*conf.peak_dtime <= 0* ||*/ conf.peak_dmz <= 0 || 
     /*conf.min_xiclen <= 0 ||*/ conf.xic_width <= 0) {
    QMessageBox::warning(this, tr("Data Load"), tr("Invalid XIC parameters."));
    return;
  }

  // Collect grouping parameters.
  conf.skip_normalize_label_free = checkSkipNormalizeLabelFree->checkState() == Qt::Checked;
  //conf.align_dtime = editAligndtime->text().toDouble();
  conf.min_conds = editMinConds->text().toInt();
  conf.min_reps = editMinReps->text().toInt();
  conf.min_iruns = editMinIruns->text().toInt();
  conf.group_dtime = editGroupdtime->text().toDouble();

  // Validate grouping paramters.
  if(conf.min_conds <= 0 || conf.align_dtime <= 0  || conf.group_dtime < 0) {
    QMessageBox::warning(this, tr("Data Load"), tr("Invalid grouping parameters."));
    return;
  }

  // Collect memory related flags.
  conf.run_algorithms = checkRunAlgo->checkState() == Qt::Checked;
  conf.keep_raw = checkRawPeaks->checkState() == Qt::Checked;
  conf.keep_filtered = checkFilteredPeaks->checkState() == Qt::Checked;
  conf.load_ms2 = checkLoadMS2->checkState() == Qt::Checked;
  conf.keep_ms2 = checkKeepMS2->checkState() == Qt::Checked;

  // Collect isotope related parameters
  conf.max_label_count = spinMaxLabel->value();
  conf.skip_normalize_ratios = checkSkipNormalizeIsotope->checkState() == Qt::Checked;
  conf.isotope_15N = checkIsotope15N->checkState() == Qt::Checked;
  conf.isotope_by_search = checkUngroupedIsotopes->checkState() == Qt::Checked;

  // Collect MS/MS related parameters.
  conf.mstol = editMSTol->text().toDouble();
  conf.ms2tol = editMS2Tol->text().toDouble();
  if(comboMS2Units->currentIndex() == 1) conf.ms2tol_usedaltons = true;
  else conf.ms2tol_usedaltons = false;

  conf.missed_cleaves = editMissedCleaves->text().toInt();

  conf.enzyme_id = max(0, comboEnzyme->currentIndex());
  conf.ms2_fdr = editFDR->text().toDouble();
  conf.fixed_15N = checkFixed15N->checkState() == Qt::Checked;

  conf.acetyl_stoich_plus5 = checkAcetylPlus5->checkState() == Qt::Checked;

  // Depending on which items were selected, set the (fixed) modifications as active or not
  for(size_t i = 0; i < conf.fxMods.size(); i++) {
    if(listFixedMods->item(i)->isSelected()) conf.fxMods[i].active = true;
    else conf.fxMods[i].active = false;
  }
  
  // Do the same for variable modifications.
  for(size_t i = 0; i < conf.varMods.size(); i++) {
    if(listVarMods->item(i)->isSelected()) conf.varMods[i].active = true;
    else conf.varMods[i].active = false;
  }
  
  // Get active heavy isotopes.
  conf.activeIsotopes.clear();
  for(size_t idx = 0; idx < conf.isotopes.size(); idx++) 
    if(listIsotopes->item(idx)->isSelected()) conf.activeIsotopes.push_back(idx);

  // Collect MS/MS modification related parameters.
  conf.max_mod_aa = editMaxModAA->text().toInt();

  // Validate fixed modifications. Only one fixed modification per amino acid!
  /*for(size_t i = 0; i < conf.fxMods.size(); i++) {
    conf.varMods[i].active
    }*/

  // If everything is OK, accept and close dialog.
  accept();
}

void MSConfigWin::PopulateFixedMods() {
  listFixedMods->clear();
  for(size_t i = 0; i < conf.fxMods.size(); i++) {
    listFixedMods->addItem(FixedModText(conf.fxMods[i]));
    listFixedMods->item(i)->setSelected(conf.fxMods[i].active);
  }
}

QString MSConfigWin::FixedModText(fxmod_t &mod) {
  QString text = QString::number(mod.deltaMass, 'f', 10);
  text += "  " + QString(mod.AAs.c_str()) + "  " + QString(mod.description.c_str());
  return text;
}

// Adds an item to the fixed item list. Validate input paramters.
void MSConfigWin::on_addFixedButton_clicked(bool) {
  string description = editFixedDescr->text().toStdString(), AAs = editFixedAAs->text().toStdString();
  double deltaMass = editFixedShift->text().toDouble();
  if(description == "" || AAs == "" || deltaMass == 0) return;

  fxmod_t mod; // Create a new fixed modification.
  mod.description = description;  mod.AAs = AAs;  mod.deltaMass = deltaMass;
  mod.active = mod.system = false;
  conf.fxMods.push_back(mod);

  listFixedMods->addItem(FixedModText(mod));
  editFixedDescr->clear(); editFixedAAs->clear(); editFixedShift->clear();
}

// Removes an item from item list.
void MSConfigWin::on_removeFixedButton_clicked(bool) {
  int item = listFixedMods->currentRow();
  if(item < 0) return; if(conf.fxMods[item].system) return;

  delete listFixedMods->takeItem(item);  // This removes the item, but we have to delete it.

  // Shift the rest of the mods left in the modificaiton array.
  while(item < (int)conf.fxMods.size() - 1) { conf.fxMods[item] = conf.fxMods[item+1]; item++; }
  conf.fxMods.pop_back();
}

QString MSConfigWin::IsotopeText(isotope_t &iso) {
  QString text = QString::number(iso.shift_heavy, 'f', 10);
  text += "  " + QString::number(iso.shift_medium, 'f', 10);
  text += "  " + QString(iso.AAs.c_str());
  text += "  " + QString(iso.description.c_str());
  return text;
}

// Fills isotope table with availible isotopes.
void MSConfigWin::PopulateIsotopes() {
  listIsotopes->clear(); 
  for(size_t i = 0; i < conf.isotopes.size(); i++) listIsotopes->addItem(IsotopeText(conf.isotopes[i]));
  for(size_t i = 0; i < conf.activeIsotopes.size(); i++) {
    if(0 <= conf.activeIsotopes[i] && conf.activeIsotopes[i] < (int)conf.isotopes.size()) {
      listIsotopes->item(conf.activeIsotopes[i])->setSelected(true);
    }
  }
}

// Adds an item to the fixed item list.
void MSConfigWin::on_addIsotopeButton_clicked(bool) {
  string description = editIsotopeDescr->text().toStdString();
  string AAs = editIsotopeAAs->text().toStdString();
  QString strShiftHeavy = editIsoShiftHeavy->text(), strShiftMedium = editIsoShiftMedium->text();
  if(description == "" || strShiftHeavy == "" || strShiftMedium == "") return;
  double shiftHeavy = strShiftHeavy.toDouble(), shiftMedium = strShiftMedium.toDouble();
  if( (shiftHeavy == 0 && shiftMedium == 0) || shiftHeavy == 0) return;
  
  isotope_t iso;
  iso.description = description; iso.AAs = AAs;
  iso.shift_heavy = shiftHeavy;  iso.shift_medium = shiftMedium;
  iso.system = false;

  conf.isotopes.push_back(iso);
  conf.updateIsotopeIdx();
  listIsotopes->addItem(IsotopeText(iso));  
  editIsotopeDescr->clear(); editIsotopeAAs->clear(); editIsoShiftHeavy->clear(); editIsoShiftMedium->clear();
}

// Removes an item from item list.
void MSConfigWin::on_removeIsotopeButton_clicked(bool) {
  int item = listIsotopes->currentRow();
  if(item < 0) return; 
  if(conf.isotopes[item].system) return;
  int row = listIsotopes->currentRow();
  delete listIsotopes->takeItem(item);  // This removes the item, but we have to delete it.
  if(listIsotopes->count() > 0) 
    listIsotopes->scrollToItem(listIsotopes->item(max(0,row-1)));

  // Shift the rest of the mods left in the modificaiton array.
  while(item < (int)conf.isotopes.size() - 1) { conf.isotopes[item] = conf.isotopes[item+1]; item++; }
  conf.isotopes.pop_back();
  conf.updateIsotopeIdx();
}


// Fill list display with variable modifications.
void MSConfigWin::PopulateVarMods() {
  listVarMods->clear();
  for(size_t i = 0; i < conf.varMods.size(); i++) {
    varmod_t &vr = conf.varMods[i];
    listVarMods->addItem(VarModText(vr));
    listVarMods->item(i)->setSelected(vr.active);
  }
}

QString MSConfigWin::VarModText(varmod_t &vr) {
  QString text;
  text += QString::number(vr.d_prec_mass, 'f', 10) + " " + QString::number(vr.max_count);
  text += " " + QString(vr.description.c_str());
  text += " " + QString(vr.abbrev.c_str());
  text += "\n";
  for(size_t i = 0; i < vr.fragmods.size(); i++) {
    text += " [" + QString(vr.fragmods[i].abbreviation.c_str()) + ",";
    text +=  QString(vr.fragmods[i].AAs.c_str()) + ",";
    text +=  QString::number(vr.fragmods[i].dfrag_mass, 'f', 10) + "]";
  }
  return text;
}

void MSConfigWin::on_addFragShiftButton_clicked(bool) {
  QString AAs = editFragAAs->text(), strFragShift = editFragShift->text(), abbrev = editAbbrev->text();
  if(AAs == "" || strFragShift == "") return;
  // Replace any pipes with a space.
  for(int i = 0; i < abbrev.length(); i++) if(abbrev[i] == '|') abbrev[i] = ' ';
  // Then take entries and add as strings to list.
  listFragShifts->addItem(abbrev + "|" + strFragShift + "|" + AAs);
  editFragAAs->setText(""); 
  editFragShift->setText(""); 
  editAbbrev->setText("");
}

void MSConfigWin::on_removeFragShiftButton_clicked(bool) {
  int item = listFragShifts->currentRow();
  if(item < 0) return; 
  delete listFragShifts->takeItem(item);  // This removes the item, but we have to delete it.
}

// Adds a variable modification if the user entered everything correctly
void MSConfigWin::on_addVarButton_clicked(bool) {
  string description = editVarDescr->text().trimmed().toStdString();
  string abbrev = editVarAbbrev->text().trimmed().toStdString();
  QString strPMass = editDeltaPMass->text();

  // Handle cases where nothing useful was entered.
  if(listFragShifts->count() == 0) return;
  if(description == "" || strPMass == ""  || abbrev == "") return;

  double d_prec_mass = strPMass.toDouble(); if(d_prec_mass == 0) return;

  varmod_t vr;  
  vr.description = description;
  vr.abbrev = abbrev;
  vr.d_prec_mass = d_prec_mass;
  vr.max_count = spinMaxCount->value();
  vr.active = false;  vr.system = false;

  for(int i = 0; i < listFragShifts->count(); i++) {
    QString shift_AAs = listFragShifts->item(i)->text();
    QString abbrev = shift_AAs.section("|", 0, 0);
    QString shift = shift_AAs.section("|", 1, 1);
    QString AAs = shift_AAs.section("|", 2, 2);
    fragmod_t fm;
    fm.abbreviation = abbrev.toStdString();
    // Note that I calculate the dfrag mass by considering the
    // precursor mass here.
    fm.dfrag_mass = d_prec_mass + shift.toDouble();
    fm.AAs = AAs.toStdString();
    vr.fragmods.push_back(fm);
  }

  conf.varMods.push_back(vr);
  conf.updateVarModIdx();
  listVarMods->addItem(VarModText(vr));  

  editVarDescr->setText(""); editDeltaPMass->setText(""); listFragShifts->clear();
}

void MSConfigWin::on_removeVarButton_clicked(bool) {
  int item = listVarMods->currentRow();
  if(item < 0) return; if(conf.varMods[item].system) return;

  int row = listVarMods->currentRow();
  delete listVarMods->takeItem(item);  // This removes the item, but we have to delete it.  
  if(listVarMods->count() > 0) 
    listVarMods->scrollToItem(listIsotopes->item(max(0,row-1)));


  // Shift the rest of the mods left in the modificaiton array.
  while(item < (int)conf.varMods.size() - 1) { conf.varMods[item] = conf.varMods[item+1]; item++; }

  conf.varMods.pop_back();
  conf.updateVarModIdx();
}

