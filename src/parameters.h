#ifndef READ_PARAMETERS_H
#define READ_PARAMETERS_H
 
#include <string>
#include <vector>
#include <cassert>

#include <tclap/CmdLine.h>

using namespace std;

class Parameters {
protected:
    int    dim;
    string mesh_path;
    string field_path;
    int    field_type;
    bool   FGLOBAL;
    string fieldG_path;
    int    n_regions;
    int    isoval_type;
    vector<double> isoval_vals;
    bool   DENOISE;
    bool   ISOCONTOURS;
    bool   ANALYZE;
    bool   CLEAN;
    bool   SMOOTH;
    int    n_iter;
    double clean_thresh;
    bool   SIGMA;
    string out_path;
    int    out_level;
    bool   GUI;
    bool   VERBOSE;

public:
    Parameters(){}
    ~Parameters(){}

    void read_cmdline( TCLAP::CmdLine &cmd, int argc, char *argv[] ) {
        TCLAP::ValueArg<string> mesh_path_arg   ("m","mesh_path","path of the input mesh.",true,"","string",cmd);
        TCLAP::ValueArg<string> field_path_arg  ("f","field_path","path of the input field.",true,"","string",cmd);
        TCLAP::ValueArg<string> fieldG_path_arg ("g","fieldG_path","path of the input global field.",false,"","string",cmd);
        TCLAP::ValueArg<string> out_path_arg    ("o","out_path","path to a directory where to save all outputs.",false,"","string",cmd);

        TCLAP::ValueArg<int> dim_arg         ("d","dim","dimension of the problem (2 or 3).",true,2,"integer",cmd);
        TCLAP::ValueArg<int> field_type_arg  ("t","field_type","input field defined on cells (1) or vertices (2).",true,1,"integer",cmd);
        TCLAP::ValueArg<int> n_regions_arg   ("r","n_regions","number of regions to be computed in the domain.",true,1,"integer",cmd);
        TCLAP::ValueArg<int> isoval_type_arg ("i","isoval_type","(1) equispaced, (2) percentiles, (3) assigned.",false,1,"integer",cmd);
        TCLAP::ValueArg<int> n_iter_arg      ("n","n_iter","max number of iterations.",false,1,"integer",cmd);
        TCLAP::ValueArg<int> out_level_arg   ("l","out_level","output level: (1) domain, (2) regions, (3) subregions.",false,0,"integer",cmd);

        TCLAP::ValueArg<double> clean_thresh_arg ("e","clean_thresh","percentual min size of the regions.",false,0,"double",cmd);
        TCLAP::MultiArg<double> isoval_vals_arg ("v","isoval_vals","isovalues (percentiles or explicit).",false,"vector<double>",cmd);

        TCLAP::SwitchArg FGLOBAL_switch     ("G","FGLOBAL","use of the global field.",cmd,false);
        TCLAP::SwitchArg DENOISE_switch     ("D","DENOISE","use the denoised field for isoregions.",cmd,false);
        TCLAP::SwitchArg ISOCONTOURS_switch ("I","ISOCONTOURS","compute isocontours/isosurfaces.",cmd,false);
        TCLAP::SwitchArg ANALYZE_switch     ("A","ANALYZE","analyze regions.",cmd,false);
        TCLAP::SwitchArg CLEAN_switch       ("C","CLEAN","cleaning of small regions.",cmd,false);
        TCLAP::SwitchArg SMOOTH_switch      ("S","SMOOTH","smoothing of the boundaries.",cmd,false);
        TCLAP::SwitchArg SIGMA_switch       ("M","SIGMA","use standard deviation.",cmd,false);
        TCLAP::SwitchArg GUI_switch         ("U","GUI","launch graphical interface.",cmd,false);
        TCLAP::SwitchArg VERBOSE_switch     ("V","verbose","print debug information.",cmd,false);

        cmd.parse( argc, argv );

        if (mesh_path_arg.isSet())    mesh_path = mesh_path_arg;
        if (field_path_arg.isSet())   field_path = field_path_arg;
        if (fieldG_path_arg.isSet())  fieldG_path = fieldG_path_arg;
        if (out_path_arg.isSet())     out_path = out_path_arg;
        if (out_level_arg.isSet())    out_level = out_level_arg;

        if (dim_arg.isSet())          dim = dim_arg.getValue();
        if (field_type_arg.isSet())   field_type = field_type_arg.getValue();
        if (n_regions_arg.isSet())    n_regions = n_regions_arg.getValue();
        if (isoval_type_arg.isSet())  isoval_type = isoval_type_arg.getValue();
        if (n_iter_arg.isSet())       n_iter = n_iter_arg.getValue();
        if (n_iter_arg.isSet())       n_iter = n_iter_arg.getValue();

        if (clean_thresh_arg.isSet()) clean_thresh = clean_thresh_arg.getValue();
        if (isoval_vals_arg.isSet())  isoval_vals = isoval_vals_arg.getValue();

        FGLOBAL = FGLOBAL_switch.getValue();
        DENOISE = DENOISE_switch.getValue();
        ISOCONTOURS = ISOCONTOURS_switch.getValue();
        ANALYZE = ANALYZE_switch.getValue();
        CLEAN = CLEAN_switch.getValue();
        SMOOTH = SMOOTH_switch.getValue();
        SIGMA = SIGMA_switch.getValue();
        GUI = GUI_switch.getValue();
        VERBOSE = VERBOSE_switch.getValue();

        assert((dim==2 || dim==3) && "Error: unsupported dim");
        assert((field_type==1 || field_type==2) && "Error: unsupported field_type");
        assert((isoval_type==1 || isoval_type==2 || isoval_type==3) && "Error: unsupported isoval_type");
        if (isoval_type > 1)
            assert((int)isoval_vals.size()==n_regions+1 && "Error: wrong number of isoval_vals");
        assert((out_level==1 || out_level==2 || out_level==3) && "Error: unsupported out_level");
    }

    void print_pars() {
        cout << "begin  print_pars"   << endl;
        cout << "dim: "               << dim << endl;
        cout << "mesh_path: "         << mesh_path << endl;
        cout << "field_path: "        << field_path << endl;
        cout << "field_type: "        << field_type << endl;
        cout << "FGLOBAL: "           << FGLOBAL << endl;
        cout << "fieldG_path: "       << fieldG_path << endl;
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
        cout << "out_path: "          << out_path << endl;
        cout << "out_level: "         << out_level << endl;
        cout << "GUI: "               << GUI << endl;
        cout << "VERBOSE: "           << VERBOSE << endl;
        cout << "end of print_pars\n"   << endl;
    }

    string get_MESH_PATH ()       { return mesh_path; }
    string get_FIELD_PATH ()      { return field_path; }
    int    get_FIELD_TYPE ()      { return field_type; }
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
    int    get_OUT_LEVEL ()       { return out_level; }
    bool   get_GUI ()             { return GUI; }
    bool   get_ISOCONTOURS ()     { return ISOCONTOURS; }
    bool   get_VERBOSE ()         { return VERBOSE; }
};

#endif // READ_PARAMETERS_H
