#ifndef READ_PARAMETERS_H
#define READ_PARAMETERS_H

/* class for reading the input parameters contained in the file 'parameters.run' */
 
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>

#include <tclap/CmdLine.h>

using namespace std;

class Parameters {
private:
  static istream & eatline(istream & s) {
    while ( s.get() != '\n' && s.good() ) {}
    return s ;
  }

protected:
    int    dim;             // dimension of the problem
    string mesh_path;       // path of the input mesh
    string field_path;      // path of the input field
    int    field_type;      // input field defined on cells (1) or vertices (2)
    string field2_path;     // path of the input field #2
    bool   FGLOBAL;         // use of the global field
    string fieldG_path;     // path of the input global field
    string out_path;        // path to a directory where to save all outputs
    int    n_regions;       // number of regions to be computed in the domain
    int    isoval_type;     // (1) equispaced, (2) percentiles, (3) assigned
    vector<double> isoval_vals; // isovalues (percentiles or explicit)
    bool   DENOISE;         // use the denoised field for isoregions
    bool   ISOCONTOURS;     // compute isocontours/isosurfaces
    bool   ANALYZE;         // analyze regions
    bool   CLEAN;           // cleaning of small regions
    bool   SMOOTH;          // smoothing of the boundaries
    int    n_iter;          // max number of iterations
    double clean_thresh;    // percentual min size of the regions
    bool   SIGMA;           // use standard deviation
    bool   gui;             // launch graphical interface
    bool   verbose;         // print debug information

    void set_defaults() {
        mesh_path     = "";
        field_path    = "";
        field_type    = 1;
        field2_path   = "";
        fieldG_path   = "";
        FGLOBAL       = false;
        DENOISE       = false;
        dim           = 2;
        n_regions     = 1;
        isoval_type   = 1;
        isoval_vals   = {};
        ANALYZE       = false;
        CLEAN         = false;
        SMOOTH        = false;
        n_iter        = 0;
        clean_thresh  = 1.;
        SIGMA         = false;
        out_path      = ".";
        gui           = true;
        ISOCONTOURS   = false;
        verbose       = false;
    }

    void read_mesh_path(ifstream & inpf)    { string s;
                                              inpf >> s;
                                              mesh_path = string(HOME_PATH) + s; }
    void read_field_path(ifstream & inpf)   { string s;
                                              inpf >> s;
                                              field_path = string(HOME_PATH) + s; }
    void read_field_type(ifstream & inpf)   { inpf >> field_type; }
    void read_field2_path(ifstream & inpf)  { string s;
                                              inpf >> s;
                                              field2_path = string(HOME_PATH) + s; }
    void read_fieldG_path(ifstream & inpf)  { string s;
                                              inpf >> s;
                                              fieldG_path = string(HOME_PATH) + s; }
    void read_FGLOBAL(ifstream & inpf)      { inpf >> FGLOBAL; }
    void read_DENOISE(ifstream & inpf)      { inpf >> DENOISE; }
    void read_dim(ifstream & inpf)          { inpf >> dim; }
    void read_n_regions(ifstream & inpf)    { inpf >> n_regions; }
    void read_isoval_type(ifstream & inpf)  { inpf >> isoval_type; }
    void read_isoval_vals(ifstream & inpf)  { isoval_vals.clear();
                                              string token;
                                              while (inpf >> token) {
                                                  if (token == "#") break;
                                                  try { isoval_vals.push_back(stod(token)); }
                                                  catch (const invalid_argument& e) { cerr << "Invalid input: " << token << endl; }
                                              }
                                              if (isoval_type > 1)
                                                  assert((int)isoval_vals.size()==n_regions+1 && "Error: wrong number of isoval_vals");
                                            }
    void read_ANALYZE(ifstream & inpf)      { inpf >> ANALYZE; }
    void read_CLEAN(ifstream & inpf)        { inpf >> CLEAN; }
    void read_SMOOTH(ifstream & inpf)       { inpf >> SMOOTH; }
    void read_n_iter(ifstream & inpf)       { inpf >> n_iter; }
    void read_clean_thresh(ifstream & inpf) { inpf >> clean_thresh; }
    void read_SIGMA(ifstream & inpf)        { inpf >> SIGMA; }
    void read_out_path(ifstream & inpf)     { inpf >> out_path; }
    void read_gui(ifstream & inpf)          { inpf >> gui; }
    void read_ISOCONTOURS(ifstream & inpf)  { inpf >> ISOCONTOURS; }
    void read_verbose(ifstream & inpf)      { inpf >> verbose; }

    bool read_line( string keywd, ifstream & inpf ) {
      bool retval = true;
      if      ( keywd == "mesh_path" )      { read_mesh_path(inpf); }
      else if ( keywd == "field_path" )     { read_field_path(inpf); }
      else if ( keywd == "field_type" )     { read_field_type(inpf); }
      else if ( keywd == "field2_path" )    { read_field2_path(inpf); }
      else if ( keywd == "fieldG_path" )    { read_fieldG_path(inpf); }
      else if ( keywd == "FGLOBAL" )        { read_FGLOBAL(inpf); }
      else if ( keywd == "DENOISE" )        { read_DENOISE(inpf); }
      else if ( keywd == "dim" )            { read_dim(inpf); }
      else if ( keywd == "n_regions" )      { read_n_regions(inpf); }
      else if ( keywd == "isoval_type" )    { read_isoval_type(inpf); }
      else if ( keywd == "isoval_vals" )    { read_isoval_vals(inpf); }
      else if ( keywd == "ANALYZE" )        { read_ANALYZE(inpf); }
      else if ( keywd == "CLEAN" )          { read_CLEAN(inpf); }
      else if ( keywd == "SMOOTH" )         { read_SMOOTH(inpf); }
      else if ( keywd == "n_iter" )         { read_n_iter(inpf); }
      else if ( keywd == "clean_thresh" )   { read_clean_thresh(inpf); }
      else if ( keywd == "SIGMA" )          { read_SIGMA(inpf); }
      else if ( keywd == "out_path" )       { read_out_path(inpf); }
      else if ( keywd == "gui" )            { read_gui(inpf); }
      else if ( keywd == "ISOCONTOURS" )    { read_ISOCONTOURS(inpf); }
      else if ( keywd == "verbose" )        { read_verbose(inpf); }
      else                                  { retval = false ; }
      if ( retval )                         { inpf >> eatline; }
      return retval;
    }

public:
    Parameters(){ set_defaults(); }
    ~Parameters(){}

    void set_MESH_PATH      (const string value)  { mesh_path = value; }
    void set_FIELD_PATH     (const string value)  { field_path = value; }
    void set_FIELD_TYPE     (const int value)     { field_type = value; }
    void set_FIELD2_PATH    (const string value)  { field2_path = value; }
    void set_FIELDG_PATH    (const string value)  { fieldG_path = value; }
    void set_FGLOBAL        (const bool value)    { FGLOBAL = value; }
    void set_DENOISE        (const bool value)    { DENOISE = value; }
    void set_DIM            (const int value)     { dim = value; }
    void set_N_REGIONS      (const int value)     { n_regions = value; }
    void set_ISOVAL_TYPE    (const int value)     { isoval_type = value; }
    void set_ISOVAL_VALS    (const vector<double> value) { isoval_vals = value; }
    void set_ANALYZE        (const bool value)    { ANALYZE = value; }
    void set_CLEAN          (const bool value)    { CLEAN = value; }
    void set_SMOOTH         (const bool value)    { SMOOTH = value; }
    void set_N_ITER         (const int value)     { n_iter = value; }
    void set_CLEAN_THRESH   (const double value)  { clean_thresh = value; }
    void set_SIGMA          (const bool value)    { SIGMA = value; }
    void set_OUT_PATH       (const string value)  { out_path = value; }
    void set_GUI            (const bool value)    { gui = value; }
    void set_ISOCONTOURS    (const bool value)    { ISOCONTOURS = value; }
    void set_VERBOSE        (const bool value)    { verbose = value; }

    string get_MESH_PATH ()       { return mesh_path; }
    string get_FIELD_PATH ()      { return field_path; }
    int    get_FIELD_TYPE ()      { return field_type; }
    string get_FIELD2_PATH ()     { return field2_path; }
    string get_FIELDG_PATH ()     { return fieldG_path; }
    bool   get_FGLOBAL ()         { return FGLOBAL; }
    bool   get_DENOISE ()         { return DENOISE; }
    int    get_DIM ()             { return dim; }
    int    get_N_REGIONS ()       { return n_regions; }
    int    get_ISOVAL_TYPE ()     { return isoval_type; }
    vector<double> get_ISOVAL_VALS() { return isoval_vals; }
    bool   get_ANALYZE ()         { return ANALYZE; }
    bool   get_CLEAN ()           { return CLEAN; }
    bool   get_SMOOTH ()          { return SMOOTH; }
    int    get_N_ITER ()          { return n_iter; }
    double get_CLEAN_THRESH ()    { return clean_thresh; }
    bool   get_SIGMA ()           { return SIGMA; }
    string get_OUT_PATH ()        { return out_path; }
    bool   get_GUI ()             { return gui; }
    bool   get_ISOCONTOURS ()     { return ISOCONTOURS; }
    bool   get_VERBOSE ()         { return verbose; }

    void print_pars() {
      cout << "begin  print_pars"   << endl;
      cout << "dim: "               << dim << endl;
      cout << "mesh_path: "         << mesh_path << endl;
      cout << "field_path: "        << field_path << endl;
      cout << "field_type: "        << field_type << endl;
      cout << "field2_path: "       << field2_path << endl;
      cout << "FGLOBAL: "           << FGLOBAL << endl;
      cout << "fieldG_path: "       << fieldG_path << endl;
      cout << "out_path: "          << out_path << endl;
      cout << "n_regions: "         << n_regions << endl;
      cout << "isoval_type: "       << isoval_type << endl;
      for (uint i=0; i<isoval_vals.size(); ++i)
          cout << "isoval_vals(" << i << "): " << isoval_vals.at(i) << endl;
      cout << "DENOISE: "           << DENOISE << endl;
      cout << "ISOCONTOURS: "       << ISOCONTOURS << endl;
      cout << "ANALYZE: "           << ANALYZE << endl;
      cout << "CLEAN: "             << CLEAN << endl;
      cout << "SMOOTH: "            << SMOOTH << endl;
      cout << "n_iter: "            << n_iter << endl;
      cout << "clean_thresh: "      << clean_thresh << endl;
      cout << "SIGMA: "             << SIGMA << endl;
      cout << "gui: "               << gui << endl;
      cout << "verbose: "           << verbose << endl;
      cout << "end of print_pars\n"   << endl;
    }
};


class Read_Parameters : public Parameters {

private:
  static istream & get_comments( istream & is ) {
    bool newline = true;
    while( newline ) {
      newline = false;
      while( !is.eof() && (is.peek()==' ' || is.peek()=='\t' || is.peek()=='\n') ) { is.get(); }
      while( !is.eof() && is.peek()=='#' ) {
        while( !is.eof() && is.peek()!='\n' ) { is.get(); }
        newline = true;
      }
    }
    return is;
  }

  string fname;
  void fatal_error( string fname ) {
    cerr << "fatal error in opening file " << fname << endl << flush;
    exit(0);
  }

  void read_line( ifstream & inpf ) {
    string keywd = "";
    inpf >> get_comments >> keywd;
    Parameters::read_line( keywd, inpf );
  }

public:
  Read_Parameters() {}
  Read_Parameters( string _fname ) : fname(_fname) {}
  ~Read_Parameters() {}

  void print_parameters() {
    Parameters::print_pars();
  }

  void read_file( string _fname="" ) {
    string inp_fname = _fname=="" ? fname : _fname;
    ifstream inpf( inp_fname.c_str() );
    if ( !inpf.good() )
        fatal_error(inp_fname);
    while ( !inpf.eof() )
        read_line( inpf );
    inpf.close() ;
  }

  void read_cmdline( TCLAP::CmdLine &cmd, int argc, char *argv[] ) {
      TCLAP::ValueArg<string> mesh_path_arg   ("m","mesh_path","path of the input mesh.",true,"","string", cmd);
      TCLAP::ValueArg<string> field_path_arg  ("f","field_path","path of the input field.",true,"","string", cmd);
      TCLAP::ValueArg<string> fieldG_path_arg ("g","fieldG_path","path of the input global field.",false,"","string", cmd);
      TCLAP::ValueArg<string> out_path_arg    ("o","out_path","path to a directory where to save all outputs.",false,"","string", cmd);

      TCLAP::ValueArg<int> dim_arg         ("d","dim","dimension of the problem (2 or 3).",true,2,"integer", cmd);
      TCLAP::ValueArg<int> field_type_arg  ("t","field_type","input field defined on cells (1) or vertices (2).",true,1,"integer", cmd);
      TCLAP::ValueArg<int> n_regions_arg   ("r","n_regions","number of regions to be computed in the domain.",true,1,"integer", cmd);
      TCLAP::ValueArg<int> isoval_type_arg ("i","isoval_type","(1) equispaced, (2) percentiles, (3) assigned.",false,1,"integer", cmd);
      TCLAP::ValueArg<int> n_iter_arg      ("n","n_iter","max number of iterations.",false,1,"integer", cmd);

      TCLAP::ValueArg<double> clean_thresh_arg ("e","clean_thresh","percentual min size of the regions.",false,0,"double", cmd);
      TCLAP::MultiArg<double> isoval_vals_arg ("v","isoval_vals","isovalues (percentiles or explicit).",false,"vector<double>", cmd);

      TCLAP::SwitchArg FGLOBAL_switch     ("G","FGLOBAL","use of the global field.", cmd, false);
      TCLAP::SwitchArg DENOISE_switch     ("D","DENOISE","use the denoised field for isoregions.", cmd, false);
      TCLAP::SwitchArg ISOCONTOURS_switch ("I","ISOCONTOURS","compute isocontours/isosurfaces.", cmd, false);
      TCLAP::SwitchArg ANALYZE_switch     ("A","ANALYZE","analyze regions.", cmd, false);
      TCLAP::SwitchArg CLEAN_switch       ("C","CLEAN","cleaning of small regions.", cmd, false);
      TCLAP::SwitchArg SMOOTH_switch      ("S","SMOOTH","smoothing of the boundaries.", cmd, false);
      TCLAP::SwitchArg SIGMA_switch       ("M","SIGMA","use standard deviation.", cmd, false);
      TCLAP::SwitchArg gui_switch         ("U","gui","launch graphical interface.", cmd, false);
      TCLAP::SwitchArg verbose_switch     ("V","verbose","print debug information.", cmd, false);

      cmd.parse( argc, argv );

      if (mesh_path_arg.isSet())    mesh_path = mesh_path_arg;
      if (field_path_arg.isSet())   field_path = field_path_arg;
      if (fieldG_path_arg.isSet())  fieldG_path = fieldG_path_arg;
      if (out_path_arg.isSet())     out_path = out_path_arg;

      if (dim_arg.isSet())          dim = dim_arg.getValue();
      if (field_type_arg.isSet())   field_type = field_type_arg.getValue();
      if (n_regions_arg.isSet())    n_regions = n_regions_arg.getValue();
      if (isoval_type_arg.isSet())  isoval_type = isoval_type_arg.getValue();
      if (n_iter_arg.isSet())       n_iter = n_iter_arg.getValue();

      if (clean_thresh_arg.isSet()) clean_thresh = clean_thresh_arg.getValue();
      if (isoval_vals_arg.isSet())  isoval_vals = isoval_vals_arg.getValue();

      if (FGLOBAL_switch.isSet())   FGLOBAL = FGLOBAL_switch.getValue();
      if (DENOISE_switch.isSet())   DENOISE = DENOISE_switch.getValue();
      if (ISOCONTOURS_switch.isSet()) ISOCONTOURS = ISOCONTOURS_switch.getValue();
      if (ANALYZE_switch.isSet())   ANALYZE = ANALYZE_switch.getValue();
      if (CLEAN_switch.isSet())     CLEAN = CLEAN_switch.getValue();
      if (SMOOTH_switch.isSet())    SMOOTH = SMOOTH_switch.getValue();
      if (SIGMA_switch.isSet())     SIGMA = SIGMA_switch.getValue();
      if (gui_switch.isSet())       gui = gui_switch.getValue();
      if (verbose_switch.isSet())   verbose = verbose_switch.getValue();

      assert((dim==2 || dim==3) && "Error: unsupported dim");
      assert((field_type==1 || field_type==2) && "Error: unsupported field_type");
      assert((isoval_type==1 || isoval_type==2 || isoval_type==3) && "Error: unsupported isoval_type");
      if (isoval_type > 1)
          assert((int)isoval_vals.size()==n_regions+1 && "Error: wrong number of isoval_vals");
  }
};

#endif // READ_PARAMETERS_H
