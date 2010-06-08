# Here are the known FORTRAN name mangling schemes:
#   gfortran:  __module_name__function_name
#              __module_name_MOD_function_name
#   portland   module_name_function_name_
#   intel      module_name_mp_function_name_
#   pathscale  TEST_ROUTINE.in.TEST_MODULE
#   cray       test_routine$test_module_
#
# Probably there are as many others as there are FORTRAN compilers. :(
#

AC_DEFUN([__FC_MODULE_NAME_MANGLING],
[_AC_FORTRAN_ASSERT()dnl
AC_CACHE_CHECK([for Fortran module name-mangling scheme],
               ac_cv_[]_AC_LANG_ABBREV[]_module_mangling,
[AC_COMPILE_IFELSE(
[      module test_module
      contains
      subroutine foo_bar()
      return
      end subroutine foo_bar
      end module test_module],
[mv conftest.$ac_objext cfortran_test.$ac_objext

  ac_save_LIBS=$LIBS
  LIBS="cfortran_test.$ac_objext $LIBS $[]_AC_LANG_PREFIX[]LIBS"

  AC_LANG_PUSH(C)dnl
  ac_success=no
  for ac_foobar in gnu:__test_module__foo_bar gnu44:__test_module_MOD_foo_bar portland:test_module_foo_bar_ intel:test_module_mp_foo_bar_ pathscale:TEST_ROUTINE.in.FOO_BAR cray:foo_bar$test_module_ ; do
      ac_func="`echo $ac_foobar | cut -d: -f2`"
      AC_LINK_IFELSE([AC_LANG_CALL([], [$ac_func])],
                     [ac_success=yes; break 1])
  done
  AC_LANG_POP(C)dnl
  ac_cv_[]_AC_LANG_ABBREV[]_module_mangling=`echo $ac_foobar | cut -d: -f1`
  LIBS=$ac_save_LIBS
  rm -f cfortran_test* conftest*],
  [AC_MSG_FAILURE([cannot compile a simple Fortran program])])
])
])# __FC_MODULE_NAME_MANGLING

# _FC_MODULE_NAME_MANGLING
# ----------------------
AC_DEFUN([_FC_MODULE_NAME_MANGLING],
[AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])dnl
AC_REQUIRE([AC_FC_DUMMY_MAIN])dnl
AC_LANG_PUSH(Fortran)dnl
__FC_MODULE_NAME_MANGLING
AC_LANG_POP(Fortran)dnl
])# _FC_MODULE_NAME_MANGLING


# _FC_MODULE_WRAPPERS
# ---------------
# Defines C macros {F77,FC}_MOD_FUNC(name,NAME) and {F77,FC}_MOD_FUNC_(name,NAME) to
# properly mangle the names of C identifiers, and C identifiers with
# underscores, respectively, so that they match the name mangling
# scheme used by the Fortran compiler.
AC_DEFUN([_FC_MODULE_WRAPPERS],
[_AC_FORTRAN_ASSERT()dnl
AH_TEMPLATE(_AC_FC[_MOD_FUNC],
    [Define to a macro mangling the given C identifier (in lower and upper
     case), which must not contain underscores, for linking with Fortran.])dnl
AH_TEMPLATE(_AC_FC[_MOD_FUNC_],
    [As ]_AC_FC[_MOD_FUNC, but for C identifiers containing underscores.])dnl
case $ac_cv_[]_AC_LANG_ABBREV[]_module_mangling in
  "gnu")
          AC_DEFINE(_AC_FC[_MOD_FUNC(mod,MOD,name,NAME)],  [__##mod##__##name])
          AC_DEFINE(_AC_FC[_MOD_FUNC_(mod,MOD,name,NAME)],  [__##mod##__##name]) ;;
  "gnu44")
          AC_DEFINE(_AC_FC[_MOD_FUNC(mod,MOD,name,NAME)],  [__##mod##_MOD_##name])
          AC_DEFINE(_AC_FC[_MOD_FUNC_(mod,MOD,name,NAME)],  [__##mod##_MOD_##name]) ;;
  "portland")
          AC_DEFINE(_AC_FC[_MOD_FUNC(mod,MOD,name,NAME)],  [mod##_##name##_])
          AC_DEFINE(_AC_FC[_MOD_FUNC_(mod,MOD,name,NAME)],  [mod##_##name##_]) ;;
  "intel")
          AC_DEFINE(_AC_FC[_MOD_FUNC(mod,MOD,name,NAME)],  [mod##_mp_##name##_])
          AC_DEFINE(_AC_FC[_MOD_FUNC_(mod,MOD,name,NAME)],  [mod##_mp_##name##_]) ;;
  "pathscale")
          AC_DEFINE(_AC_FC[_MOD_FUNC(mod,MOD,name,NAME)],  [NAME##.in.##MOD])
          AC_DEFINE(_AC_FC[_MOD_FUNC_(mod,MOD,name,NAME)],  [NAME##.in.##MOD]) ;;
  "cray")
          AC_DEFINE(_AC_FC[_MOD_FUNC(mod,MOD,name,NAME)],  [name##$##mod##_])
          AC_DEFINE(_AC_FC[_MOD_FUNC_(mod,MOD,name,NAME)],  [name##$##mod##_]) ;;
  *)
          AC_MSG_WARN([unknown Fortran name-mangling scheme])
          ;;
esac
])# _FC_MODULE_WRAPPERS

# FC_MODULE_WRAPPERS
# --------------
AC_DEFUN([FC_MODULE_WRAPPERS],
[AC_REQUIRE([_FC_MODULE_NAME_MANGLING])dnl
AC_LANG_PUSH(Fortran)dnl
_FC_MODULE_WRAPPERS
AC_LANG_POP(Fortran)dnl
])# FC_MODULE_WRAPPERS
