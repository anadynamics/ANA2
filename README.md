# Analysis of Null Areas 

        ANA input.pdb -c config_file.cfg -d trajectory.nc -o output.pdb 


Precise cavity definition. Flexible tracking. Fast performance.
Accurate volume calculation. Ugly output.

![](/abstract_fig_cut.png)
![](./cpp_logo.png)

## Installation

See installation info on the manual.

## Usage example

Check the quickstart example.


## Release History

* 1.0.0
    * Initial release.

* 1.0.1
    * Volume output can be redirected to a file.

* Road to 2.0.0
    * Major refactoring.
    * Switched to C++ 17.
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
