# SlicerDRRGenerator

A 3D Slicer extension that generates Digitally Reconstructed Radiographs (DRRs).

Note that there are alternative methods/modules for DRR computation:
- DRRs can be created using the `Volume Rendering` Slicer core module, by using the `CT-X-Ray` preset adn setting 3D view background to black. This method is useful for real-time rendering of DRR images, but it is limited to computing 8 bit per pixel output.
- `DRR Image Computation` module in SlicerRT extension can compute DRR image with the imaging geometry specified by a radiation therapy beam.
