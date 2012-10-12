#!/usr/bin/env perl

use strict;
use warnings;

use Math::GSL::SF qw/gsl_sf_bessel_J1/;

use constant {
  epsilon_0 => 8.85e-12,
  mu_0      => 1.26e-6,
  qe        => 1.4e-19,
  pi        => 4*atan2(1,1),
  bess_j_zero_01 => 2.405,
};

use constant { steve => gsl_sf_bessel_J1(bess_j_zero_01)**2 };

use PDL;
use PDL::FFT qw/fft/;
use PDL::Graphics::Prima::Simple;

my $res_freq = 2.8838e9;
my $dielec_const = 10.4;
my $pulse_duration = 4.0e-12;
my $speed = 1.0e8;
my $pulse_width_rho = 1e-3;
my $electrons = 2000;
my $loop_area = 1.26e-5;
my $q_factor = 6000;

my $charge = qe * $electrons;
my $permittivity = $dielec_const * epsilon_0;
my $pulse_width_z = $pulse_duration * $speed;
my $pulse_volume = pi * $pulse_width_z * $pulse_width_rho**2;
my $charge_density = $charge / $pulse_volume;
my $res_freq_rad = 2 * pi * $res_freq;

sub max_length {
  my $res = shift;
  return $speed / ( 2 * $res );
}

sub max_radius {
  my $res = shift;
  my $denom = 2 * pi * $res * sqrt( $permittivity * mu_0 );
  return bess_j_zero_01 / $denom;
}

sub flight_time {
  my $res = shift;
  return max_length( $res ) / $speed;
}

sub Emax {
  my $res = shift;
  my $l = max_length( $res );
  my $r = max_radius( $res );
  my $sqrt = sqrt( ( pi * $l * $r )**2 * $permittivity * steve );
  return 1 / $sqrt;
}

sub pulse_train {
  my ($z, $t, $duration) = @_;
  my $arg = ($z - $speed * $t) / ( $duration * $speed / 2 );
  return $charge_density * $speed * exp( - ($arg ** 2) );
}

my $tmax = 1e-9;
my $steps = 1e3;
my $t = $tmax * (xvals( 2*$steps + 1 ) - $steps) / $steps;

my $train = pulse_train( 0, $t, 4e-12 );
my $imag = $train->zeros;
fft($train, $imag);
line_plot $t, abs($train**2 + $imag**2);

