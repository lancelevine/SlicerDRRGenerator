cmake_minimum_required(VERSION 3.5)

project(DRRGenerator)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://slicerdrrregistration.000webhostapp.com/drrgenerator")
set(EXTENSION_CATEGORY "Filtering")
set(EXTENSION_CONTRIBUTORS "Lance Levine")
set(EXTENSION_DESCRIPTION "This extension allows the user to generate DRRs via a highly customizable GUI")
set(EXTENSION_ICONURL "https://slicerdrrregistration.000webhostapp.com/wp-content/uploads/2020/01/DRRGenerator.png")
set(EXTENSION_SCREENSHOTURLS "https://slicerdrrregistration.000webhostapp.com/wp-content/uploads/2020/01/screenshot-1024x556.png")
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
