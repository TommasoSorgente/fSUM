#ifndef READ_PARAMETERS_H
#define READ_PARAMETERS_H

/* class for reading the input parameters contained in the file 'parameters.run' */
 
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>

using namespace std;

class Parameters {
private:
  static istream & eatline(istream & s) {
    while ( s.get() != '\n' && s.good() ) {}
    return s ;
  }

protected:
    string mesh_path;       // path of the input mesh
    string field1_path;     // path of the input field #1
    string field2_path;     // path of the input field #2
    string fieldG_path;     // path of the input field #2
    int    dim;             // dimension of the problem
    bool   smooth_field;    // smooth the field values by averaging
    int    n_regions;       // number of regions to be computed in the domain
    vector<double> isovals; // isovalues
    bool   cut_mesh;        // cut mesh along isovalues
    bool   log_mode;        // use the logarithm of the field values (for non-normalized data)
    int    n_iter;          // max number of iterations
    double filter_thresh;   // min area of the regions (percentual)
    int    ref_type;        // 1 for edge split, 2 for cell split
    double shp_thresh;      // refine regions with sqrt(area)/perimeter smaller than this
    double dim_thresh;      // refine regions with area smaller than this (percentual)
    double var_thresh;      // refine regions with quality variance smaller than this
    string out_path;        // path do a directory where to save all outputs
    string out_format;      // file format for saving the output mesh
    bool   gui;             // launch graphical interface
    bool   show_iso;        // show isocontours/isosurfaces on the final plot
    bool   verbose;         // print debug information

    void set_defaults() {
        mesh_path     = "";
        field1_path   = "";
        field2_path   = "";
        fieldG_path   = "";
        dim           = 2;
        smooth_field  = false;
        n_regions     = 1;
        isovals       = {};
        cut_mesh      = false;
        log_mode      = false;
        n_iter        = 0;
        filter_thresh = 1.;
        ref_type      = 2;
        shp_thresh    = 1.;
        dim_thresh    = 1.;
        var_thresh    = 1.;
        out_path      = ".";
        out_format    = "";
        gui           = true;
        show_iso      = false;
        verbose       = false;
    }

    void read_mesh_path(ifstream & inpf)    { string s;
                                              inpf >> s;
                                              mesh_path = string(HOME_PATH) + s; }
    void read_field1_path(ifstream & inpf)  { string s;
                                              inpf >> s;
                                              field1_path = string(HOME_PATH) + s; }
    void read_field2_path(ifstream & inpf)  { string s;
                                              inpf >> s;
                                              field2_path = string(HOME_PATH) + s; }
    void read_fieldG_path(ifstream & inpf)  { string s;
                                              inpf >> s;
                                              fieldG_path = string(HOME_PATH) + s; }
    void read_dim(ifstream & inpf)          { inpf >> dim; }
    void read_smooth_field(ifstream & inpf) { inpf >> smooth_field; }
    void read_n_regions(ifstream & inpf)    { inpf >> n_regions; }
    void read_isovals(ifstream & inpf)      { isovals.clear();
                                              while ( inpf.get()!='\n' && inpf.good() ) {
                                                double f;
                                                inpf >> f;
                                                isovals.push_back(f);
                                              }
                                              assert((int)isovals.size()==n_regions+1 && "Error: wrong number of isovals"); }
    void read_cut_mesh(ifstream & inpf)     { inpf >> cut_mesh; }
    void read_log_mode(ifstream & inpf)     { inpf >> log_mode; }
    void read_n_iter(ifstream & inpf)       { inpf >> n_iter; }
    void read_filter_thresh(ifstream & inpf){ inpf >> filter_thresh; }
    void read_ref_type(ifstream & inpf)     { inpf >> ref_type; }
    void read_shp_thresh(ifstream & inpf)   { inpf >> shp_thresh; }
    void read_dim_thresh(ifstream & inpf)   { inpf >> dim_thresh; }
    void read_var_thresh(ifstream & inpf)   { inpf >> var_thresh; }
    void read_out_path(ifstream & inpf)     { inpf >> out_path; }
    void read_out_format(ifstream & inpf)   { inpf >> out_format; }
    void read_gui(ifstream & inpf)          { inpf >> gui; }
    void read_show_iso(ifstream & inpf)     { inpf >> show_iso; }
    void read_verbose(ifstream & inpf)      { inpf >> verbose; }

    bool read_line( string keywd, ifstream & inpf ) {
      bool retval = true;
      if      ( keywd == "mesh_path" )      { read_mesh_path(inpf); }
      else if ( keywd == "field1_path" )    { read_field1_path(inpf); }
      else if ( keywd == "field2_path" )    { read_field2_path(inpf); }
      else if ( keywd == "fieldG_path" )    { read_fieldG_path(inpf); }
      else if ( keywd == "dim" )            { read_dim(inpf); }
      else if ( keywd == "smooth_field" )   { read_smooth_field(inpf); }
      else if ( keywd == "n_regions" )      { read_n_regions(inpf); }
      else if ( keywd == "isovals" )        { read_isovals(inpf); }
      else if ( keywd == "cut_mesh" )       { read_cut_mesh(inpf); }
      else if ( keywd == "log_mode" )       { read_log_mode(inpf); }
      else if ( keywd == "n_iter" )         { read_n_iter(inpf); }
      else if ( keywd == "filter_thresh" )  { read_filter_thresh(inpf); }
      else if ( keywd == "ref_type" )       { read_ref_type(inpf); }
      else if ( keywd == "shp_thresh" )     { read_shp_thresh(inpf); }
      else if ( keywd == "dim_thresh" )     { read_dim_thresh(inpf); }
      else if ( keywd == "var_thresh" )     { read_var_thresh(inpf); }
      else if ( keywd == "out_path" )       { read_out_path(inpf); }
      else if ( keywd == "out_format" )     { read_out_format(inpf); }
      else if ( keywd == "gui" )            { read_gui(inpf); }
      else if ( keywd == "show_iso" )       { read_show_iso(inpf); }
      else if ( keywd == "verbose" )        { read_verbose(inpf); }
      else                                  { retval = false ; }
      if ( retval )                         { inpf >> eatline; }
      return retval;
    }

public:
    Parameters(){ set_defaults(); }
    ~Parameters(){}

    void set_MESH_PATH      (const string value)  { mesh_path = value; }
    void set_FIELD1_PATH    (const string value)  { field1_path = value; }
    void set_FIELD2_PATH    (const string value)  { field2_path = value; }
    void set_FIELDG_PATH    (const string value)  { fieldG_path = value; }
    void set_DIM            (const int value)     { dim = value; }
    void set_SMOOTH_FIELD   (const bool value)    { smooth_field = value; }
    void set_N_REGIONS      (const int value)     { n_regions = value; }
    void set_ISOVALS (const vector<double> value) { isovals = value; }
    void set_CUT_MESH       (const bool value)    { cut_mesh = value; }
    void set_LOG_MODE       (const bool value)    { log_mode = value; }
    void set_N_ITER         (const int value)     { n_iter = value; }
    void set_FILTER_THRESH  (const double value)  { filter_thresh = value; }
    void set_REF_TYPE       (const int value)     { ref_type = value; }
    void set_SHP_THRESH     (const double value)  { shp_thresh = value; }
    void set_DIM_THRESH     (const double value)  { dim_thresh = value; }
    void set_VAR_THRESH     (const double value)  { var_thresh = value; }
    void set_OUT_PATH       (const string value)  { out_path = value; }
    void set_OUT_FORMAT     (const string value)  { out_format = value; }
    void set_GUI            (const bool value)    { gui = value; }
    void set_SHOW_ISO       (const bool value)    { show_iso = value; }
    void set_VERBOSE        (const bool value)    { verbose = value; }

    string get_MESH_PATH ()       { return mesh_path; }
    string get_FIELD1_PATH ()     { return field1_path; }
    string get_FIELD2_PATH ()     { return field2_path; }
    string get_FIELDG_PATH ()     { return fieldG_path; }
    int    get_DIM ()             { return dim; }
    bool   get_SMOOTH_FIELD ()    { return smooth_field; }
    int    get_N_REGIONS ()       { return n_regions; }
    vector<double> get_ISOVALS () { return isovals; }
    bool   get_CUT_MESH ()        { return cut_mesh; }
    bool   get_LOG_MODE ()        { return log_mode; }
    int    get_N_ITER ()          { return n_iter; }
    double get_FILTER_THRESH ()   { return filter_thresh; }
    int    get_REF_TYPE ()        { return ref_type; }
    double get_SHP_THRESH ()      { return shp_thresh; }
    double get_DIM_THRESH ()      { return dim_thresh; }
    double get_VAR_THRESH ()      { return var_thresh; }
    string get_OUT_PATH ()        { return out_path; }
    string get_OUT_FORMAT ()      { return out_format; }
    bool   get_GUI ()             { return gui; }
    bool   get_SHOW_ISO ()        { return show_iso; }
    bool   get_VERBOSE ()         { return verbose; }

    void print_pars() {
      cout << "begin  print_pars"   << endl;
      cout << "mesh_path: "         << mesh_path << endl;
      cout << "field1_path: "       << field1_path << endl;
      cout << "field2_path: "       << field2_path << endl;
      cout << "fieldG_path: "       << fieldG_path << endl;
      cout << "dim: "               << dim << endl;
      cout << "smooth_field: "      << smooth_field << endl;
      cout << "n_regions: "         << n_regions << endl;
      for (uint i=0; i<isovals.size(); ++i)
          cout << "isovals(" << i << ": " << isovals.at(i) << endl;
      cout << "cut_mesh: "          << cut_mesh << endl;
      cout << "log_mode: "          << log_mode << endl;
      cout << "n_iter: "            << n_iter << endl;
      cout << "filter_thresh: "     << filter_thresh << endl;
      cout << "ref_type: "          << ref_type << endl;
      cout << "shp_thresh: "        << shp_thresh << endl;
      cout << "dim_thresh: "        << dim_thresh << endl;
      cout << "var_thresh: "        << var_thresh << endl;
      cout << "out_path: "          << out_path << endl;
      cout << "out_format: "        << out_format << endl;
      cout << "gui: "               << gui << endl;
      cout << "show_iso: "          << show_iso << endl;
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
};

#endif // READ_PARAMETERS_H
