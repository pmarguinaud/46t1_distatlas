#!/usr/bin/perl -w

use strict;
use FileHandle;
use File::Path;
use File::Basename;


sub label
{
  my $grid = shift;

  if ($grid =~ m/^L(\d+)x(\d+)/o)
    {
      return "Global lat/lon 0 .. 360 ${1}x${2}";
    }

  if ($grid =~ m/^X(\d+)x(\d+)/o)
    {
      return "Global lat/lon -180 .. +180 ${1}x${2}";
    }

  if ($grid =~ m/^S(\d+)x(\d+)/o)
    {
      return "Global shifted lat/lon ${1}x${2}";
    }

  if ($grid =~ m/^N(\d+)$/o)
    {
      return "Gaussian T" . (2 * $1);
    }

  if ($grid =~ m/^fort.4.t(\d+)$/o)
    {
      return "Gaussian T$1";
    }

  if ($grid eq 'ICMSHARPEINIT')
    {
      return 'ARPEGE T32 c2.4';
    }

  if ($grid eq 'ICMSHAROMINIT')
    {
      return 'AROME 2.5km 64x64';
    }

  if ($grid =~ m/^PFFCST(\w+)\+0000$/o)
    {
      return $1;
    }

  if ($grid =~ m/^fort\.4\.(\d+)x(\d+)$/o)
    {
      return "AROME France 50km ${1}x${2}";
    }

  if ($grid eq 'fort.4.32x32_100km')
    {
      return 'AROME France 100km 32x32';
    }

  return $grid;
}


sub nn
{
  return $ENV{SLURM_NNODES} if (exists $ENV{SLURM_NNODES});
  return 1;
}

sub mpiauto
{
  my ($nnp, $openmp, @exec) = @_;

  my $mpiauto = "/home/gmap/mrpm/marguina/SAVE/mpiauto/mpiauto";

  my @mpiauto = ($mpiauto, qw (--wrap --wrap-stdeo), -nnp => $nnp, 
                 -openmp => $openmp, -nn => &nn (), '--', @exec);

  printf ("@mpiauto\n");

  local $ENV{DR_HOOK} = 0;
  local $ENV{DR_HOOK_NOT_MPI} = 1;
  die if (system ("ulimit -s unlimited; @mpiauto"));

  unlink ('linux_bind.txt');

}


my $PACK = '/home/gmap/mrpm/marguina/pack/46t1_distatlas.01.I185274INTELMPI184274MT.x';


sub lfitools
{
  local $ENV{DR_HOOK} = 0;
  local $ENV{DR_HOOK_NOT_MPI} = 1;
  my @exec = ("$PACK/bin/lfitools", @_);
  system (@exec) && die ("@exec failed\n");
}

sub catfa
{

  my $prefix = shift;

  my @f = <$prefix.fa.*>;

  if (@f)
    {
      &lfitools (qw (lfi_alt_index --lfi-file-in), @f, '--lfi-file-out', "$prefix.fa");
      &lfitools (qw (lfi_alt_pack --lfi-file-in), "$prefix.fa", '--lfi-file-out', "$prefix.pack.fa");
    }

}


sub run
{
  my ($nnp, $openmp, $option, @t) = @_;
  
  for my $t (@t)
    {
      my ($grid1, $grid2) = @$t;
      my $dir = "$grid1-$grid2";

      &rmtree ($dir);
      &mkpath ($dir);
      chdir ($dir);

      for my $f (<$PACK/data/*>)
        {
          symlink ($f, &basename ($f));
        }

      &mpiauto 
        (
          $nnp, $openmp,
          "$PACK/bin/ATLAS_ARPEGE_F", '--grid1' => $grid1,  
          '--grid2' => $grid2, '--write2', $option,
          '--light1', '--light2'
        );

      for my $p (qw (XYZ1 XYZ2 XYZ2I4 XYZ2IA))
        {
          &catfa ($p);
        }


      chdir ('..');

    }

}


my @si4 =
(
  [qw (L80x40          fort.4.t32          )], # Global lat/lon to various domains
  [qw (L80x40          ICMSHARPEINIT       )],
  [qw (L80x40          fort.4.64x64        )],
  [qw (L80x40          ICMSHAROMINIT       )],
  [qw (S80x40          fort.4.t32          )],
  [qw (L80x40          L60x30              )],
  [qw (L400x200        PFFCSTEUROC25+0000  )],
  [qw (N32             N64                 )], # ARPEGE to ARPEGE
  [qw (fort.4.t32      fort.4.64x64        )], # ARPEGE to AROME
  [qw (fort.4.64x64    fort.4.32x32        )], # AROME to AROME
);


if (0){
&mkpath ('si4');
chdir ('si4');
&run (4, 2, '--interp4', @si4);
chdir ('..');
}


my @siA =
(
  [qw (L320x160       fort.4.32x32         )], # Global lat/lon to various domains
  [qw (X160x80        fort.4.t32           )],
  [qw (L160x80        fort.4.t32           )],
  [qw (L400x200       L40x20               )],
  [qw (N80            fort.4.32x32         )], # Global Gaussian to various domains
  [qw (N320           fort.4.64x64         )],
  [qw (N320           N80                  )],
  [qw (fort.4.64x64   fort.4.32x32_100km   )], # AROME to AROME
);

if (0){
&mkpath ('siA');
chdir ('siA');
&run (4, 2, '--interpA', @siA);
chdir ('..');
}


my @bi4 =
(
  [qw (N1024          L4000x2000           )],
  [qw (N2000          L400x200             )],
  [qw (N2000          fort.4.512x512       )],
);

if (0){
&mkpath ('bi4');
chdir ('bi4');
&run (8, 16, '--interp4', @bi4);
chdir ('..');
}


my @biA =
(
  [qw (L40000x20000   N1024                )],
  [qw (L4000x2000     N512                 )],
  [qw (L40000x20000   fort.4.512x512       )],
);

if (1){
&mkpath ('biA');
chdir ('biA');
&run (8, 16, '--interpA', @biA);
chdir ('..');
}






