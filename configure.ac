## autoconf configure.ac && autoreconf -if
AC_PREREQ([2.69])
AC_INIT([glm], [0.3-0])
AC_CONFIG_SRCDIR([src/glm_fit.c])

R_SCMD="${R_HOME}/bin/Rscript -e"


# SEXPtools
SEXPTOOLS_CPPFLAGS=`${R_SCMD} "SEXPtools:::cppflags()"`
SEXPTOOLS_LDFLAGS=`${R_SCMD} "SEXPtools:::ldflags()"`


# Report
echo " "
echo "****************** Results of glm package configure *****************"
echo " "
echo ">> SEXPTOOLS_CPPFLAGS = ${SEXPTOOLS_CPPFLAGS}"
echo ">> SEXPTOOLS_LDFLAGS = ${SEXPTOOLS_LDFLAGS}"
echo " "
echo "*********************************************************************"
echo " "


# Substitute templates
AC_SUBST(SEXPTOOLS_CPPFLAGS)
AC_SUBST(SEXPTOOLS_LDFLAGS)

AC_OUTPUT(src/Makevars)

