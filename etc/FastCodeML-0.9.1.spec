# This is an RPM spec file for a FastCodeML
# 
#

# Base directory
%define category	Phylogeny
%define _prefix         /software/%{category}/%{name}/%{version}
%define module_path     /software/module/%{category}/%{name}
# Those three because RedHat defines it differently
%define _defaultdocdir %{_datadir}/doc
%define _mandir %{_datadir}/man
%define _infodir %{_datadir}/info

# Those definitions make sure the RPM build machinery also uses our binutils tools
%global __ld                   /software/bin/ld
%global __nm                   /software/bin/nm
%global __objcopy              /software/bin/objcopy
%global __objdump              /software/bin/objdump
%global __ranlib               /software/bin/ranlib
%global __strip                /software/bin/strip

BuildRoot:              %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
Summary:                CodeML's Branch Site Model
License:                GPL
Vendor:                 SIB/Vital-IT
Name:                   FastCodeML
Version:                0.9.1
Release:                1
Source:                 %{name}-%{version}.tar.gz
Prefix:                 %{_prefix}
Packager:               hstockin
Group:		        %{category}
# Require our build environment
BuildRequires: Vcluster-gcc, Vcluster-gcc-c++
BuildRequires: Vcluster-libstdc++-devel
BuildRequires: Vcluster-OSlinks
BuildRequires: Vcluster-gcc-gfortran
BuildRequires: Vcluster-libgfortran
BuildRequires: Vcluster-libgfortran-static
BuildRequires: nlopt
BuildRequires: OpenBLAS
BuildRequires: Vcluster-cmake
BuildRequires: Vcluster-bzip2
BuildRequires: Vcluster-bzip2-libs
BuildRequires: Vcluster-bzip2-devel
BuildRequires: Vcluster-boost
BuildRequires: Vcluster-boost-build
BuildRequires: Vcluster-boost-chrono
BuildRequires: Vcluster-boost-date-time
BuildRequires: Vcluster-boost-devel
BuildRequires: Vcluster-boost-filesystem
BuildRequires: Vcluster-boost-graph
BuildRequires: Vcluster-boost-iostreams
BuildRequires: Vcluster-boost-jam
BuildRequires: Vcluster-boost-math
BuildRequires: Vcluster-boost-program-options
BuildRequires: Vcluster-boost-regex
BuildRequires: Vcluster-boost-serialization
BuildRequires: Vcluster-boost-signals
BuildRequires: Vcluster-boost-static
BuildRequires: Vcluster-boost-system
BuildRequires: Vcluster-boost-test
BuildRequires: Vcluster-boost-thread

%description

%prep
export PATH=/software/bin${PATH:+:$PATH}
%setup -q


%build
# Put the PATH to our tools first
export PATH=/software/bin${PATH:+:$PATH}

# Assume Vital-IT installation pathes
export BLAS_LIB_DIR="/software/Utility/OpenBLAS/0.2.5/lib"
export LAPACK_LIB_DIR="/software/Utility/OpenBLAS/0.2.5/include"
export NLOPT_LIB_DIR="/software/Utility/nlopt/2.3/lib"
export NLOPT_INCLUDE_DIR="/software/Utility/nlopt/2.3/include"
export MATH_LIB_NAMES="openblas;lapack;gfortranbegin;gfortran"
export CXX=g++

cmake .
make

# Create the module file parametrically
cat << \EOF > %{name}-%{version}.module
#%Module1.0
##
## %{name} modulefile
##
proc ModulesHelp { } {
        global version

        puts stderr "\tThis module loads the %{name} environment"
        puts stderr "\tVersion $version"
        puts stderr "\n\tType 'module list' to list all the loaded modules"
        puts stderr "\tand 'module avail' to list all the availables ones."
}

module-whatis   "%{name} version %{version}"

# for Tcl script use only
set     version     %{version}
set     root        %{_prefix}

prepend-path    PATH              	$root/bin
EOF

%check
# Put the PATH to our tools first
export PATH=/software/bin${PATH:+:$PATH}


%install
# Put the PATH to our tools first
export PATH=/software/bin${PATH:+:$PATH}


rm -rf ${RPM_BUILD_ROOT}
mkdir -p ${RPM_BUILD_ROOT}%{_prefix}/bin
install fast ${RPM_BUILD_ROOT}%{_prefix}/bin

# install module file
mkdir -p ${RPM_BUILD_ROOT}%{module_path}
install %{name}-%{version}.module ${RPM_BUILD_ROOT}%{module_path}/%{version}


%clean
rm -rf $RPM_BUILD_ROOT


%files
%defattr(-,root,root,-)
%{_prefix}/bin/*

%{module_path}/%{version}


%changelog
* Thu Feb 14 2013 Heinz Stockinger <Heinz.Stockinger@isb-sib.ch> - 0.9.1
Created.
