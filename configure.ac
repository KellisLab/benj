AC_INIT([benj], [0.0.7])

AC_MSG_NOTICE([Creating VERSION file])
echo m4_defn([AC_PACKAGE_VERSION]) > VERSION

AC_CHECK_PROG([R_FOUND], [R], [yes], [no])
if test "x$R_FOUND" = xno; then
  AC_MSG_ERROR([R is required but it's not installed.])
fi
AC_CHECK_PROG([WGET_FOUND], [wget], [yes], [no])
if test "x$WGET_FOUND" = xno; then
  AC_MSG_ERROR([wget is required but it's not installed.])
fi
AC_PROG_MAKE_SET

# Automake initialization
AM_INIT_AUTOMAKE([foreign])

# Define the directory where files should be downloaded
TENX_DIR="inst/extdata/10x"
AC_SUBST(TENX_DIR)

# Create the extdata directory if it does not exist
AC_MSG_CHECKING([for inst/extdata directory])
AS_MKDIR_P([${TENX_DIR}])
AC_MSG_RESULT([created])

# Download files using wget if they don't exist
AC_MSG_CHECKING([for required files])
REQUIRED_FILES="Dual_Index_Kit_TT_Set_A.csv Dual_Index_Kit_NN_Set_A.csv Dual_Index_Kit_NT_Set_A.csv Single_Index_Kit_T_Set_A.csv Single_Index_Kit_N_Set_A.csv"
for file in $REQUIRED_FILES; do
  if test ! -f "${TENX_DIR}/$file"; then
    AC_MSG_RESULT([missing $file, downloading...])
    wget -P ${TENX_DIR} --no-check-certificate "https://cdn.10xgenomics.com/raw/upload/v1655151897/support/in-line documents/$file"
  else
    AC_MSG_RESULT([found $file])
  fi
done

AC_CANONICAL_HOST
# Detect OS and substitute the appropriate compilers
AC_MSG_CHECKING([for operating system $host_os])
case "$host_os" in
  *darwin*)
    AC_MSG_RESULT([Mac OS X detected])
    AC_SUBST([CC], [conda-forge::clang_osx-64])
    AC_SUBST([CXX], [conda-forge::clangxx_osx-64])
    AC_SUBST([FC], [conda-forge::gfortran_osx-64])
    ;;
  *linux*)
    AC_MSG_RESULT([Linux detected])
    AC_SUBST([CC], [conda-forge::gcc_linux-64])
    AC_SUBST([CXX], [conda-forge::gxx_linux-64])
    AC_SUBST([FC], [conda-forge::gfortran_linux-64])
    ;;
  *)
    AC_MSG_ERROR([unsupported operating system])
    ;;
esac

AC_CONFIG_FILES([Makefile conda/benj.yml DESCRIPTION])
AC_OUTPUT
