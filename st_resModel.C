void st_resModel( bool redo=true, TString model="generic" ){
   gSystem->CompileMacro("myDSCBPdf.cxx","k");
   gSystem->CompileMacro("resModel.C","k");
   resModel( redo, model );
}
