dnl -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(UWPKFCVM, 1.1.0, software@scec.org)
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([foreign tar-pax no-exeext])

AC_PROG_MKDIR_P
LT_INIT
AC_PROG_INSTALL
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MKDIR_P

dnl TAR
AC_PATH_PROG(TAR, tar)
if test -z "$TAR" ; then
  AC_MSG_FAILURE([cannot find 'tar' program.])
  TAR=`echo "Error: tar is not installed." ; false`
fi

dnl PROJ
AC_ARG_ENABLE([proj],
    [AC_HELP_STRING([--enable-proj],
        [install Proj @<:@default=yes@:>@])],
        [if test "$enableval" = "yes" ; then install_proj=yes; else install_proj=no; fi],
        [install_proj=yes])
AM_CONDITIONAL([INSTALL_PROJ], [test "$install_proj" = yes])
AC_ARG_WITH([proj-incdir],
    [AC_HELP_STRING([--with-proj-incdir], [directory containing Proj header files @<:@default=no@:>@])],
        [with_proj_incdir="$withval"],
        [with_proj_incdir=no])
AC_SUBST(with_proj_incdir)
AC_ARG_WITH([proj-libdir],
    [AC_HELP_STRING([--with-proj-libdir], [directory containing Proj library @<:@default=no@:>@])],
        [with_proj_libdir="$withval"],
        [with_proj_libdir=no])
AC_SUBST(with_proj_libdir)
dnl
dnl Setup environment so dependencies are used in build
dnl


CFLAGS="$CFLAGS"
CPPFLAGS="-I$prefix/include $CPPFLAGS"
LDFLAGS="-L$prefix/lib -L$prefix/lib64 $LDFLAGS"

dnl
dnl Verify configuration
dnl

##check optional large data path
##UWPKFCVM_LARGEDATA_DIR=$CVM_LARGEDATA_DIR/model/uwpkfcvm
if test x"$CVM_LARGEDATA_DIR" != x; then
   # test directory existence
   UWPKFCVM_LARGEDATA_DIR=$CVM_LARGEDATA_DIR/model/uwpkfcvm
# In docker container building.. this is nto accessible  yet
   if test x"$CVM_IN_DOCKER" != x; then
     AM_CONDITIONAL([WITH_UWPKFCVM_LARGEDATA_DIR], true)
     AC_SUBST(UWPKFCVM_LARGEDATA_DIR)
   else
     AC_CHECK_FILE([$UWPKFCVM_LARGEDATA_DIR/vp.dat], 
                 [AM_CONDITIONAL([WITH_UWPKFCVM_LARGEDATA_DIR], true) AC_SUBST(UWPKFCVM_LARGEDATA_DIR)],
                 [AM_CONDITIONAL(WITH_UWPKFCVM_LARGEDATA_DIR, false)])
   fi
else
   AM_CONDITIONAL(WITH_UWPKFCVM_LARGEDATA_DIR, false)
fi

AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  data/Makefile
  test/Makefile
  ])

AC_OUTPUT

AM_CONDITIONAL([BUILD_SHARED_LIBRARY], [test "$enable_shared" = yes])
if test "$enable_shared" = yes; then
  AC_MSG_NOTICE([Building shared libraries. UCVM must be built with shared libraries for plugins to work.])
else
  AC_MSG_WARN([Building static libraries. Shared libraries are not enabled.])
fi


dnl End of file
