# Analysis of Null Areas 

        ANA input.pdb -c config_file.cfg -d trajectory.nc -f output.pdb 


Precise cavity definition. Flexible tracking. Fast performance.
Accurate volume calculation. Ugly output.

![](/abstract_fig_cut.png)
![](./cpp_logo.png)

## Installation

https://anadynamics.netlify.com/docs/install_instructions.html

## Usage example

https://anadynamics.netlify.com/docs/quickstart.html

## Release History

* 1.0.0
    * Initial release.

* 1.0.1
    * Volume output can be redirected to a file.
    * ANA now appends ".pdb" to output PDB files if the user didn't do it already.

* Road to 2.0.0
    * Major refactoring.
    * Switched to C++ 17.
    * New on Static:
        * New output_wall terminal option to set output filename.
    * New on NDD:
        * No need to input displaced PDBs. ANA will read Amber modes or a raw
text file with vectors and move the inpud pdb atoms forward and backwards along
those vectors.
        * ANA will output either, the volumes of the perturbed cavities, the
numerical derivative, or the flexibility coefficient for that cavity.
        * Default scaling_factors are better adjusted to prevent cavity changes from
high frequency modes from being overestimated.
    * Deprecated:
        * Cell filtering by facet area.
    * Deprecations for NDD:
        * No wall residues/atoms output.
        * included area is now mandatory.
        * Fixed to high precision mode.
        * No cavity output.
    * Deprecations for MD:
    * Deprecations for Static:
        * Pymol CGO output.

### Contact info
[@gpbarletta](https://twitter.com/gpbarletta) - pbarletta@gmail.com

### License
Distributed under the [GPL](https://www.gnu.org/copyleft/gpl.html) license.
