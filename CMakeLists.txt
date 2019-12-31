cmake_minimum_required(VERSION 3.5)

project(DRRGenerator)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/DRRGenerator")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Lance Levine")
set(EXTENSION_DESCRIPTION "This extension allows the user to generate DRRs via a highly customizable GUI")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/DRRGenerator.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/DRRGenerator/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(DRRGeneratorModule)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
