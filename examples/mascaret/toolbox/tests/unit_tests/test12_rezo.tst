// This file is released under the 3-clause BSD license. See COPYING-BSD.
//=================================
toolbox_dir=getenv("toolbox_dir");
c = filesep();

// creation du modele
[erreur, id] = MASCARET_create();
assert_checkequal(id,1);

// importation du modele
path_xml = "file://"+toolbox_dir+c+".."+c+"test"+c+"Test_Plan"+c+"Test12"+c+"data"+c+"REZO"+c+"xml"+c+"mascaret0.xcas";
TabNomFichier = [
                 strsubst(path_xml,'\','/'), ..
                 toolbox_dir+c+".."+c+"test"+c+"Test_Plan"+c+"Test12"+c+"data"+c+"REZO"+c+"xml"+c+"mascaret0.geo", ..
                 toolbox_dir+c+".."+c+"test"+c+"Test_Plan"+c+"Test12"+c+"data"+c+"REZO"+c+"xml"+c+"mascaret0_0.loi", ..
                 toolbox_dir+c+".."+c+"test"+c+"Test_Plan"+c+"Test12"+c+"data"+c+"REZO"+c+"xml"+c+"mascaret0_1.loi", ..
                 toolbox_dir+c+"mascaret0.lis", ..
                 toolbox_dir+c+"mascaret0_ecr.opt"];
 
TypeNomFichier = ["xcas","geo","loi","loi","listing","res"];
impression = 0;
erreur = MASCARET_importModel(id,TabNomFichier,TypeNomFichier,impression);
assert_checkequal(erreur,0);

// initialisation
erreur = MASCARET_initStateName(id,toolbox_dir+c+".."+c+"test"+c+"Test_Plan"+c+"Test12"+c+"data"+c+"REZO"+c+"xml"+c+"mascaret0.lig",impression);
assert_checkequal(erreur,0);

// acces aux pas de temps de simulation
[erreur,pasTps] = MASCARET_getDouble(id,"Model.DT",0,0,0);
assert_checkequal(erreur,0);
[erreur,T0] = MASCARET_getDouble(id,"Model.InitTime",0,0,0);
assert_checkequal(erreur,0);
[erreur,TF] = MASCARET_getDouble(id,"Model.MaxCompTime",0,0,0);
assert_checkequal(erreur,0);
TF = 1000.;

[erreur,Zr] = MASCARET_getDouble(id,"Model.Zbot",1,0,0); // cote de fond de la 1ere section
assert_checkequal(erreur,0);
ZTP = zeros(101,1);
[erreur,ZTP(1)] = MASCARET_getDouble(id,"State.Z",1,0,0); // section amont
assert_checkequal(erreur,0);

tpsCalcul = pasTps;
i = 2;
// calcul
while (tpsCalcul <= TF)
  erreur = MASCARET_compute(id,T0,tpsCalcul,pasTps,impression);
  assert_checkequal(erreur,0);
  T0 = tpsCalcul;
  tpsCalcul = tpsCalcul + pasTps;
  [erreur,ZTP(i)] = MASCARET_getDouble(id,"State.Z",1,0,0); // section amont
  assert_checkequal(erreur,0);
  i = i + 1;
end

// recuperation des resultats
[erreur,nbSec,taille2,taille3] = MASCARET_getSizeVar(id,"Model.X", 0);
assert_checkequal(erreur,0);

// recuperation des resultats
Z = zeros(nbSec,1);
Q = zeros(nbSec,1);
for i = 1:nbSec
    [erreur,Z(i)] = MASCARET_getDouble(id,"State.Z",i,0,0);
    assert_checkequal(erreur,0);
    [erreur,Q(i)] = MASCARET_getDouble(id,"State.Q",i,0,0);
    assert_checkequal(erreur,0);
end

ResRef = read(toolbox_dir+c+".."+c+"test"+c+"Test_Plan"+c+"Test12"+c+"ref"+c+"res.txt",nbSec,7);

// test de la solution sur la cote
code_retour = assert_checkalmostequal(Z,ResRef(:,2),%eps,1.D-3);
assert_checktrue(code_retour);

// test de la solution sur le debit
code_retour = assert_checkalmostequal(Q,ResRef(:,3),%eps,1.D-3);
assert_checktrue(code_retour);

// test sur la solution analytique de la hauteur d'eau amont au temps = 200 s
code_retour = assert_checkalmostequal(ZTP(21)-Zr,4.718,%eps,1.D-3);
assert_checktrue(code_retour);

// test sur la solution analytique de la hauteur d'eau amont au temps = 400 s
code_retour = assert_checkalmostequal(ZTP(41)-Zr,4.445,%eps,2.D-3);
assert_checktrue(code_retour);

// test sur la solution analytique de la hauteur d'eau amont au temps = 600 s
code_retour = assert_checkalmostequal(ZTP(61)-Zr,4.180,%eps,3.D-3);
assert_checktrue(code_retour);

// test sur la solution analytique de la hauteur d'eau amont au temps = 800 s
code_retour = assert_checkalmostequal(ZTP(81)-Zr,3.923,%eps,3.D-3);
assert_checktrue(code_retour);

// test sur la solution analytique de la hauteur d'eau amont au temps = 1000 s
code_retour = assert_checkalmostequal(ZTP(101)-Zr,3.674,%eps,4.D-3);
assert_checktrue(code_retour);

// destruction du modele
erreur=MASCARET_delete(id);
assert_checkequal(erreur,0);

