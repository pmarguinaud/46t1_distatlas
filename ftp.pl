#!/usr/bin/perl -w
#
use strict;
use Net::FTP;

my $ftp = 'Net::FTP'->new ("meteo.fr", Debug => 1);

$ftp->login (qw (wwwguest meteo.fr));
$ftp->mkdir ("pub/marguina");
$ftp->cwd ("pub/marguina");

$ftp->binary ();

for my $f (qw (ZXYZ2.pack.fa ZXYZ2I.pack.fa ZXYZ1.pack.fa))
  {
    $ftp->put ($f, $f);
  }




