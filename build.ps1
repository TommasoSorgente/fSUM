# Get the script directory
$SCRIPTDIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$BUILDDIR = "build"
$SYSBUILDDIR = "$BUILDDIR-$env:OS"

Write-Output $SCRIPTDIR

# Navigate to the 'segmentation' directory
Set-Location -Path "$SCRIPTDIR\segmentation"
# Create the system build directory if it doesn't exist
if (-Not (Test-Path -Path $SYSBUILDDIR)) {
    New-Item -ItemType Directory -Path $SYSBUILDDIR | Out-Null
}
Set-Location -Path $SYSBUILDDIR
# Run cmake commands
cmake ..
cmake --build . --parallel 8 --config Release
Set-Location -Path ..
# Move the system build directory to the build directory
if (Test-Path -Path $BUILDDIR) {
    Remove-Item -Recurse -Force $BUILDDIR
}
Rename-Item -Path $SYSBUILDDIR -NewName $BUILDDIR

# Navigate to the 'merge_subdomains' directory
Set-Location -Path "$SCRIPTDIR\merge_subdomains"
# Create the system build directory if it doesn't exist
if (-Not (Test-Path -Path $SYSBUILDDIR)) {
    New-Item -ItemType Directory -Path $SYSBUILDDIR | Out-Null
}
Set-Location -Path $SYSBUILDDIR
# Run cmake commands
cmake ..
cmake --build . --parallel 8 --config Release
Set-Location -Path ..
# Move the system build directory to the build directory
if (Test-Path -Path $BUILDDIR) {
    Remove-Item -Recurse -Force $BUILDDIR
}
Rename-Item -Path $SYSBUILDDIR -NewName $BUILDDIR

# Return to the script directory
Set-Location -Path $SCRIPTDIR
