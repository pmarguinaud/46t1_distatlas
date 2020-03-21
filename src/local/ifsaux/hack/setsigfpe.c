#include <fenv.h>
#include <signal.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

void setsigfpe_ ();

static void handle_fpe (int signum, siginfo_t * info, void * ptr)
{

  if (info->si_code == FPE_INTDIV) fprintf (stderr, "   FPE_INTDIV     integer divide by zero          \n");
  if (info->si_code == FPE_INTOVF) fprintf (stderr, "   FPE_INTOVF     integer overflow                \n");
  if (info->si_code == FPE_FLTDIV) fprintf (stderr, "   FPE_FLTDIV     floating-point divide by zero   \n");
  if (info->si_code == FPE_FLTOVF) fprintf (stderr, "   FPE_FLTOVF     floating-point overflow         \n");
  if (info->si_code == FPE_FLTUND) fprintf (stderr, "   FPE_FLTUND     floating-point underflow        \n");
  if (info->si_code == FPE_FLTRES) fprintf (stderr, "   FPE_FLTRES     floating-point inexact result   \n");
  if (info->si_code == FPE_FLTINV) fprintf (stderr, "   FPE_FLTINV     floating-point invalid operation\n");
  if (info->si_code == FPE_FLTSUB) fprintf (stderr, "   FPE_FLTSUB     subscript out of range          \n");

  LinuxTraceBack (ptr);

  abort ();
}

void setsigfpe_ ()
{
  struct sigaction new, old;

  char * CATCHFPE = getenv ("CATCHFPE");

  if (CATCHFPE == NULL)
    return;

  feenableexcept (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

  memset (&new, 0, sizeof (new));

  new.sa_sigaction = handle_fpe;
  new.sa_flags     = SA_SIGINFO;

  sigaction (SIGFPE, &new, &old);

}

