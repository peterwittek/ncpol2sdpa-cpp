dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(ncpol2sdpa, 1.0)
AC_LANG(C++)

AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE()

CXXFLAGS="-std=c++0x"
AC_PROG_CXX

AC_COMPILE_IFELSE([AC_LANG_SOURCE(
  [[template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

    typedef check<check<bool>> right_angle_brackets;

    int a;
    decltype(a) b;

    typedef check<int> check_type;
    check_type c;
    check_type&& cr = static_cast<check_type&&>(c);]])],,
  AC_MSG_FAILURE(['$CXX $CXXFLAGS' does not accept ISO C++11]))

SYMBOLICCPP_INC=""
AC_ARG_WITH(symbolicc++-incdir,
AS_HELP_STRING([--with-symbolicc++-incdir=DIR],[SymbolicC++ include directory @<:@default /usr/include@:>@]),
[
  if test "X${withval}" = "Xno"; then :; else
    SYMBOLICCPP_INC="-I${withval}"
    AC_MSG_CHECKING(user-defined symbolicc++ includes)
    AC_MSG_RESULT([${SYMBOLICCPP_INC}])
  fi
]
)

CXXFLAGS="$CXXFLAGS $SYMBOLICCPP_INC"

SYMBOLICCPP_LIBDIR=""
AC_ARG_WITH(symbolicc++-libdir,
AS_HELP_STRING([--with-symbolicc++-libdir=DIR],[SymbolicC++ library directory @<:@default /usr/lib@:>@]),
[
  if test "X${withval}" = "Xno"; then :; else
    SYMBOLICCPP_LIBDIR="-L${withval}"
    AC_MSG_CHECKING(user-defined symbolicc++ library directory)
    AC_MSG_RESULT([${SYMBOLICCPP_LIBDIR}])

  fi
]
)
  
LIBS="$LIBS $SYMBOLICCPP_LIBDIR -lsymbolicc++"
AC_LINK_IFELSE(
      [AC_LANG_PROGRAM([#include "symbolicc++.h"],
        [Symbolic X("X")])],
      [],
      [echo "Error! You need to have symbolicc++ around."
       exit -1
      ])

AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT
