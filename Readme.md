# CRTomo - Complex Resistivity Tomography

[General description]

[a few nice inversion results]

[list of features]

[link to documentation]

[link to crtomo_tools]

## Scope

CRTomo reached a mature state and at this point new features are only added
sporadically. If you are looking for a well-tested 2D/2.5D complex resistivity
inversion code, CRTomo could be something for you. If you are new to electrical
inversion and are looking for a general purpose inversion code/framework, have
a look at PyGimli (www.pygimli.org) and its associated electrical inversion
framework BERT (https://gitlab.com/resistivity-net/bert).

## Getting help

We appreciate any suggestions, improvements, and bug reports via the github
interfaces (issues, merge requests). Due to time constraints in the in general
we cannot provide detailed usage help. However, if you have interesting
research applications and would like to use CRTomo, drop us a line and lets see
what we can do.

For technical matters please contact Maximilian Weigand
(mweigand@geo.uni-bonn.de).

For research proposals/larger collaborations, please contact Prof. A. Kemna
(kemna@geo.uni-bonn.de).


## Cite as

For scientific research we would ask you to attribute usage using the following
citation:

	Kemna, A.: Tomographic inversion of complex resistivity – theory and
	application, Ph.D.  thesis, Ruhr-Universität Bochum,
	doi:10.1111/1365-2478.12013, 2000.

For technical documentation please cite either the full git commit hash, or (if
available) a specific release with its DOI.

## License

CRTomo and its components, if not otherwise stated, is distributed under the
MIT licence. Please refer to the file **COPYING** for further information.

## Installation

If you are working with a git repository, run

	./autogen.sh

once to generate the necessary autotools files.

Otherwise, configure, compile, and install the program with

	./configure
	make
	make install

The default installation path points to $HOME/bin. You can change the
installation directory with

	./configure --prefix=$HOME/inst/

Run

	./clean_autotools.sh

to remove all auto-generated files.

## Distribution

After running ./configure, a source archive called "crtomomod*.tar.gz" can be
generated using

	make dist

## Debian Packaging

The command

	./prep_deb_package.sh

creates a Debian package in the "package/" subdirectory.

## Limiting CPU usage

CRTomoMod is compiled using openmp. The number of CPU cores that are
facilitated can be controlled using the following environment variables:

	export OMP_THREAD_LIMIT=1
	export OMP_NUM_THREADS=1
	export OPENBLAS_NUM_THREADS=1

## TODO

We should restrict the calling of external shell scripts (git version, error
messages, etc) to the ./autogen.sh script, and remove them from configure.ac.
This simplifies the Debian package generation.

