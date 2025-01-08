# Define variables

$SCRIPTDIR = Split-Path -Parent $MyInvocation.MyCommand.Path

$segmentation = "$SCRIPTDIR/build/Release/segmentation.exe"
# $segmentation = "../segmentation/build/Qt_6_8_0_for_macOS-Debug/segmentation"
# $segmentation = "../segmentation/build/Desktop_Qt_6_8_0-Debug/segmentation"

$dim = 2
$mesh_path = "$SCRIPTDIR/../data/Liguria_tri/11_SASSELLO/mesh.obj"
$field_path = "$SCRIPTDIR/../data/Liguria_tri/11_SASSELLO/field.csv"
$field_type = 1
$FGLOBAL = 0
$fieldG_path = "$SCRIPTDIR/../data/Liguria_tri/field_global.csv"

$n_regions = 5
$isoval_type = 2
$isoval_vals = @(0, 25, 50, 75, 95, 100)
$DENOISE = 1
$ISOCONTOURS = 0

$ANALYZE = 1
$CLEAN = 1
$SMOOTH = 1
$clean_thresh = 2
$n_iter = 50

$out_path = "out/"
$out_level = 1
$GUI = 1
$VERBOSE = 0

# Build the values string
$values_string = "-m $mesh_path -f $field_path -g $fieldG_path -d $dim -t $field_type -r $n_regions -i $isoval_type -n $n_iter -e $clean_thresh -o $out_path -l $out_level"

# Build the isoval string
$isoval_string = $isoval_vals -join " -v "
$isoval_string = "-v $isoval_string"


# Build the switch string
$switch_string = ""
if ($FGLOBAL -eq 1) { $switch_string += "-G " }
if ($DENOISE -eq 1) { $switch_string += "-D " }
if ($ISOCONTOURS -eq 1) { $switch_string += "-I " }
if ($ANALYZE -eq 1) { $switch_string += "-A " }
if ($CLEAN -eq 1) { $switch_string += "-C " }
if ($SMOOTH -eq 1) { $switch_string += "-S " }
if ($GUI -eq 1) { $switch_string += "-U " }
if ($VERBOSE -eq 1) { $switch_string += "-V " }

# Print the constructed command
Write-Output "launcher segmentation: $segmentation $values_string $isoval_string $switch_string"
Write-Output ""

# Run the segmentation command
Start-Process -FilePath "$segmentation" -ArgumentList "$values_string $isoval_string $switch_string" -NoNewWindow -Wait
