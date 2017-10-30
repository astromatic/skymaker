.. File Installing.rst

Installing the software
=======================

Hardware requirements
---------------------

|SkyMaker| runs in (ANSI) text-mode from a shell. A window system is
not necessary for basic operation.

The amount of memory required by |SkyMaker| depends on the size of the image
being generated. It will generally be, in bytes, of the order of 16 times the
number of pixels in the simulated image.
cores.

Obtaining |SkyMaker|
-----------------

For Linux users, the simplest way to have |SkyMaker| up and running is to install
the standard binary package the comes with your Linux distribution. Run, e.g.,
``apt-get skymaker`` (on Debian) or ``dnf install skymaker`` (Fedora) and
|SkyMaker|, as well as all its dependencies, will automatically be installed. If
you decided to install the package this way you may skip the following and move
straight to the :ref:`next section <using_skymaker>`.

However if |SkyMaker| is not available in your distribution, or to obtain the most
recent version, the |SkyMaker| source package can be downloaded from `the official
GitHub repository <https://github.com/astromatic/skymaker>`_ . One may choose
`one of the stable releases <https://github.com/astromatic/skymaker/releases>`_,
or for the fearless, `a copy of the current master development branch
<https://github.com/astromatic/skymaker/archive/master.zip>`_.

Software requirements
---------------------

|SkyMaker| has been developed on `GNU/Linux <http://en.wikipedia.org/wiki/Linux>`_
machines and should compile on any
`POSIX <http://en.wikipedia.org/wiki/POSIX>`_-compliant system (this includes
|OSX|_ and `Cygwin <http://www.cygwin.com>`_ on |Windows|_, at the price of
some difficulties with the configuration), provided that the
*development* package of the following library has been installed:

* |FFTw|_ V3.0 and above [#fftw_install]_.

On Fedora/Redhat distributions for instance, the development package above is
available as ``fftw-devel``. Note that |FFTw| is not necessary if |SkyMaker| is
linked with |Intel|'s |MKL|_ library.

Installation
------------

To install from the |GitHub| source package, you must first uncompress the
archive:

.. code-block:: console

  $ unzip skymaker-<version>.zip

A new directory called :file:`skymaker-<version>` should now appear at the current
location on your disk. Enter the directory and generate the files required by
the `autotools <http://en.wikipedia.org/wiki/GNU_Build_System>`_, which the
package relies on:

.. code-block:: console

  $ cd skymaker-<version>
  $ sh autogen.sh

A :program:`configure` script is created. This script has many options, which
may be listed with the ``--help`` option:

.. code-block:: console

  $ ./configure --help

No options are required for compiling with the default GNU C compiler
(:program:`gcc`) if all the required libraries are installed at their default
locations:

.. code-block:: console

  $ ./configure

Compared to :program:`gcc` and the librairies above, the combination of the
|Intel| compiler (:program:`icc`) and the |MKL|_ libraries can give the
|SkyMaker| executable a strong boost in performance, thanks to better
vectorized code. If :program:`icc` and the |MKL| are installed on your system
[#geticc]_ , you can take advantage of them using

.. code-block:: console

  $ ./configure --enable-mkl

Additionally, if the |SkyMaker| binary is to be run on a different machine
that does not have :program:`icc` and the |MKL| installed (e.g., a cluster
computing node), you must configure a partially statically linked executable
using

.. code-block:: console

  $ ./configure --enable-mkl --enable-auto-flags --enable-best-link

In all cases, |SkyMaker| can now be compiled with

.. code-block:: console

  $ make -j

An :file:`src/sky` executable is created. For system-wide installation, run
the usual

.. code-block:: console

  $ sudo make install

You may now check that the software is properly installed by simply
typing in your shell:

.. code-block:: console

  $ sky -v

which will return the version number and other basic information (note that
some shells require the :program:`rehash` command to be run before making a
freshly installed executable accessible in the execution path).

.. [#mac_install] Mac OS X |.dmg|_ packages should be available soon.
.. [#geticc] The Linux versions of the |Intel| compiler and |MKL| are
   `available for free to academic researchers, students, educators and open
   source contributors <http://software.intel.com/qualify-for-free-software>`_.

.. include:: keys.rst

