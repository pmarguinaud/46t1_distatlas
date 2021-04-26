#!/usr/bin/perl -w

use strict;
use FileHandle;
use File::Basename;
use File::Spec;
use Data::Dumper;
use List::Util qw (uniq);
use File::Find;

my %code;

sub slurp
{
  my $f = shift;

  unless (exists $code{$f})
    {
      $code{$f} = do { my $fh = 'FileHandle'->new ("<$f"); local $/ = undef; <$fh> };
    }

  return $code{$f};
}

sub use
{
  my $code = shift;
  my @use = ($code =~ m/^\s*use\s+(\w+)/igoms);
  return &uniq (map { lc ($_) } @use);
}

sub mod
{
  my $code = shift;
  my @mod = ($code =~ m/^\s*module\s+(\w+)/igoms);
  return &uniq (map { lc ($_) } @mod);
}

sub inc
{
  my $code = shift;
  my @inc = ($code =~ m/^\s*#include\s+"(\S+)\"/goms);
  return @inc;
}

my @src = qw (
src/local/atlas/programs/atlas-pgd-demo.F90
src/local/atlas/programs/atlas-helper.F90
src/local/atlas/interpolation/interpolation4.cc
src/local/atlas/interpolation/interpolationA.cc
src/local/atlas/fortran/interpolationA_mod.F90
src/local/atlas/fortran/interpolation4_mod.F90
src/local/atlas/fortran/gradient_mod.F90
src/local/atlas/gradient/gradient.cc
src/local/atlas/surfex/atlas-compute-aos.F90
src/local/atlas/surfex/atlas-compute-sso.F90
src/local/atlas/surfex/atlas-prep-impl.F90
src/local/atlas/surfex/atlas-pgd-impl.F90
src/local/atlas/surfex/readcovers_mod.F90
src/local/atlas/surfex/atlas-compute-covers.F90
src/local/atlas/io/atlas-fmt-null.F90
src/local/atlas/io/atlas-io-gathscat.F90
src/local/atlas/io/atlas-io-dh.F90
src/local/atlas/io/atlas-fmt.F90
src/local/atlas/io/atlas-io.F90
src/local/ifsaux/module/parkind1.F90
src/local/ifsaux/module/xrd_getoptions.F90  
src/local/ifsaux/module/xrd_unix_env.F90
src/local/ifsaux/fi_libc/fi_libc.c
src/local/ifsaux/support/iswap8.c
);

my %mod;

for my $src (grep { m/\.F90$/o } @src)
  {
    my $code = &slurp ($src);
    my @mod = &mod ($code);
    for (@mod)
      {
        $mod{$_} = $src;
      }
  }

my %inc;

&find (
  {wanted => sub { return unless ((my $f = $_) =~ m/\.h$/o); 
                   $inc{&basename ($f)} = $f; }, 
   no_chdir => 1}, 'src/local');

my $fm = 'FileHandle'->new ('>Makefile');

$fm->print (<< "EOF");

include Makefile.inc

all: pgd

EOF

my @obj;

for my $src (@src)
  {
    my $code = &slurp ($src);
    my $dir = &dirname ($src);
    $src = &basename ($src);

    my $obj;

    my @inc = uniq
              map { '-I' . 'File::Spec'->abs2rel (&dirname ($_), $dir) }
              grep { $_ } 
              map { $inc{$_} } &inc ($code);

    if ($src =~ m/\.F90$/o)
      {
        my $f90 = $src;
        ($obj = $f90) =~ s/\.F90$/.o/o;
        my @use = uniq
                  map { (my $obj = $_) =~ s/\.F90$/.o/o; $obj }
                  grep { $_ } 
                  map { $mod{$_} } 
                  &use ($code);

        push @inc, map { '-I' . 'File::Spec'->abs2rel (&dirname ($_), $dir) } @use;

        $fm->print (<< "EOF");
$dir/$obj: $dir/$src @use
	cd $dir && \$(F90) -c @inc -o $obj $f90

EOF
      }
    elsif ($src =~ m/\.cc$/o)
      {
        my $cc = $src;
        ($obj = $cc) =~ s/\.cc$/.o/o;
        $fm->print (<< "EOF");
$dir/$obj: $dir/$src 
	cd $dir && \$(CXX) -c @inc -o $obj $cc

EOF
      }
    elsif ($src =~ m/\.c$/o)
      {
        my $c = $src;
        ($obj = $c) =~ s/\.c$/.o/o;
        $fm->print (<< "EOF");
$dir/$obj: $dir/$src 
	cd $dir && \$(CC) -c @inc -o $obj $c

EOF
      }
    else
      {
        die ("Unexpected $src");
      }

    push @obj, "$dir/$obj";
  }

$fm->print (<< "EOF");

pgd: @obj
	\$(LD) -o pgd @obj \$(LIBS)

clean: 
	\\rm -f @obj pgd


EOF

$fm->close ();

