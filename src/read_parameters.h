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
    bool   mean_field;      // use the geometric mean of the field
    int    dim;             // dimension of the problem
    int    n_regions;       // number of regions to be computed in the domain
    int    isoval_type;     // (1) equispaced, (2) percentiles, (3) assigned
    vector<double> isoval_vals; // isovalues (percentiles or assigned)
    bool   cut_mesh;        // cut mesh along isovalues
    bool   filter;          // apply filtering of small regions
    bool   smoothen;        // apply smoothing of the boundaries
    int    n_iter;          // max number of iterations
    double filter_thresh;   // min area of the regions (percentual)
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
        mean_field    = false;
        dim           = 2;
        n_regions     = 1;
        isoval_type   = 1;
        isoval_vals   = {};
        cut_mesh      = false;
        filter        = false;
        smoothen      = false;
        n_iter        = 0;
        filter_thresh = 1.;
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
    void read_mean_field(ifstream & inpf)   { inpf >> mean_field; }
    void read_dim(ifstream & inpf)          { inpf >> dim; }
    void read_n_regions(ifstream & inpf)    { inpf >> n_regions; }
    void read_isoval_type(ifstream & inpf)  { inpf >> isoval_type; }
    void read_isoval_vals(ifstream & inpf)  { isoval_vals.clear();
                                              while ( inpf.get()!='\n' && inpf.good() ) {
                                                double f;
                                                inpf >> f;
                                                isoval_vals.push_back(f);
                                              }
                                              assert((int)isoval_vals.size()==n_regions+1 && "Error: wrong number of isoval_vals"); }
    void read_cut_mesh(ifstream & inpf)     { inpf >> cut_mesh; }
    void read_filter(ifstream & inpf)       { inpf >> filter; }
    void read_smoothen(ifstream & inpf)     { inpf >> smoothen; }
    void read_n_iter(ifstream & inpf)       { inpf >> n_iter; }
    void read_filter_thresh(ifstream & inpf){ inpf >> filter_thresh; }
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
      else if ( keywd == "mean_field" )     { read_mean_field(inpf); }
      else if ( keywd == "dim" )            { read_dim(inpf); }
      else if ( keywd == "n_regions" )      { read_n_regions(inpf); }
      else if ( keywd == "isoval_type" )    { read_isoval_type(inpf); }
      else if ( keywd == "isoval_vals" )    { read_isoval_vals(inpf); }
      else if ( keywd == "cut_mesh" )       { read_cut_mesh(inpf); }
      else if ( keywd == "filter" )         { read_filter(inpf); }
      else if ( keywd == "smoothen" )       { read_smoothen(inpf); }
      else if ( keywd == "n_iter" )         { read_n_iter(inpf); }
      else if ( keywd == "filter_thresh" )  { read_filter_thresh(inpf); }
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
    void set_MEAN_FIELD     (const bool value)    { mean_field = value; }
    void set_DIM            (const int value)     { dim = value; }
    void set_N_REGIONS      (const int value)     { n_regions = value; }
    void set_ISOVAL_TYPE    (const int value)     { isoval_type = value; }
    void set_ISOVAL_VALS    (const vector<double> value) { isoval_vals = value; }
    void set_CUT_MESH       (const bool value)    { cut_mesh = value; }
    void set_FILTER         (const bool value)    { filter = value; }
    void set_SMOOTHEN       (const bool value)    { smoothen = value; }
    void set_N_ITER         (const int value)     { n_iter = value; }
    void set_FILTER_THRESH  (const double value)  { filter_thresh = value; }
    void set_OUT_PATH       (const string value)  { out_path = value; }
    void set_OUT_FORMAT     (const string value)  { out_format = value; }
    void set_GUI            (const bool value)    { gui = value; }
    void set_SHOW_ISO       (const bool value)    { show_iso = value; }
    void set_VERBOSE        (const bool value)    { verbose = value; }

    string get_MESH_PATH ()       { return mesh_path; }
    string get_FIELD1_PATH ()     { return field1_path; }
    string get_FIELD2_PATH ()     { return field2_path; }
    string get_FIELDG_PATH ()     { return fieldG_path; }
    bool   get_MEAN_FIELD ()      { return mean_field; }
    int    get_DIM ()             { return dim; }
    int    get_N_REGIONS ()       { return n_regions; }
    int    get_ISOVAL_TYPE ()     { return isoval_type; }
    vector<double> get_ISOVAL_VALS () { return isoval_vals; }
    bool   get_CUT_MESH ()        { return cut_mesh; }
    bool   get_FILTER ()          { return filter; }
    bool   get_SMOOTHEN ()        { return smoothen; }
    int    get_N_ITER ()          { return n_iter; }
    double get_FILTER_THRESH ()   { return filter_thresh; }
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
      cout << "mean_field: "        << mean_field << endl;
      cout << "dim: "               << dim << endl;
      cout << "n_regions: "         << n_regions << endl;
      cout << "isoval_type: "       << isoval_type << endl;
      for (uint i=0; i<isoval_vals.size(); ++i)
          cout << "isoval_vals(" << i << "): " << isoval_vals.at(i) << endl;
      cout << "cut_mesh: "          << cut_mesh << endl;
      cout << "filter: "            << filter << endl;
      cout << "smoothen: "          << smoothen << endl;
      cout << "n_iter: "            << n_iter << endl;
      cout << "filter_thresh: "     << filter_thresh << endl;
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
