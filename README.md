# muon-angles-ee-scattering

A C script designed to run in ROOT (https://root.cern). This Monte Carlo script generates 1,000,000 events of 𝑒+𝑒−→𝜇+𝜇− scattering. 
Specifically, it uses the "Acceptance-Rejection method" to generate events that follow a specific Probability Distribution Function (PDF).

The script produces:
1. A sequence of the muon's azimuthal angles from 𝑒+𝑒− scattering that follows the PDF: F = (3/8)*[1 + (cosθ)^2].

2. A sequence of the muon's polar angles that follows the PDF of a uniform distribution.

The script creates two histograms: one for the azimuthal angle and one for the polar angle, and fits the respective PDFs. Additionally, a 2D histogram is plotted with the distributions of both the azimuthal and polar angles.

The code was developed 
**To run this script, you must have ROOT installed.**

