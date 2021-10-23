library flutter_suncalc;

import 'dart:math' as math;

const PI = math.pi;
const RAD = PI / 180;
const E = RAD * 23.4397; // obliquity of the Earth

// date/time constants and conversions
const dayMs = 1000 * 60 * 60 * 24;
const J1970 = 2440588;
const J2000 = 2451545;
const J0 = 0.0009;
final julianEpoch = DateTime.utc(-4713, 11, 24, 12, 0, 0);

class SunTiming {
  final DateTime date;
  final double lat;
  final double lng;
  DateTime? solarNoon;
  DateTime? nadir;
  DateTime? sunrise;
  DateTime? sunset;
  DateTime? sunriseEnd;
  DateTime? sunsetStart;
  DateTime? dawn;
  DateTime? dusk;
  DateTime? nauticalDawn;
  DateTime? nauticalDusk;
  DateTime? nightEnd;
  DateTime? night;
  DateTime? goldenHourEnd;
  DateTime? goldenHour;

  SunTiming(this.date, this.lat, this.lng) {
    final lw = RAD * -this.lng;
    final phi = RAD * this.lat;

    final d = toDays(this.date);
    final n = julianCycle(d, lw);
    final ds = approxTransit(0, lw, n);

    final M = solarMeanAnomaly(ds);
    final L = eclipticLongitude(M);
    final dec = declination(L, 0);

    final jnoon = solarTransitJ(ds, M, L);

    var jset, jrise;
    var diff = 0.0;

    this.solarNoon = fromJulian(jnoon);
    this.nadir = fromJulian(jnoon - 0.5);

    diff = -0.833; // 'sunrise', 'sunset'
    jset = getSetJ(diff * RAD, lw, phi, dec, n, M, L);
    jrise = jnoon - (jset - jnoon);
    this.sunrise = fromJulian(jrise);
    this.sunset = fromJulian(jset);

    diff = -0.3; //, 'sunriseEnd',    'sunsetStart' ],
    jset = getSetJ(diff * RAD, lw, phi, dec, n, M, L);
    jrise = jnoon - (jset - jnoon);
    this.sunriseEnd = fromJulian(jrise);
    this.sunsetStart = fromJulian(jset);

    diff = -6; // , 'dawn',          'dusk'        ],
    jset = getSetJ(diff * RAD, lw, phi, dec, n, M, L);
    jrise = jnoon - (jset - jnoon);
    this.dawn = fromJulian(jrise);
    ;
    this.dusk = fromJulian(jset);

    diff = -12; // , 'nauticalDawn',  'nauticalDusk'],
    jset = getSetJ(diff * RAD, lw, phi, dec, n, M, L);
    jrise = jnoon - (jset - jnoon);
    this.nauticalDawn = fromJulian(jrise);
    ;
    this.nauticalDusk = fromJulian(jset);

    diff = -18; // , 'nightEnd',      'night'       ],
    jset = getSetJ(diff * RAD, lw, phi, dec, n, M, L);
    jrise = jnoon - (jset - jnoon);
    this.nightEnd = fromJulian(jrise);
    ;
    this.night = fromJulian(jset);

    diff = 6; // , 'goldenHourEnd', 'goldenHour'  ]
    jset = getSetJ(diff * RAD, lw, phi, dec, n, M, L);
    jrise = jnoon - (jset - jnoon);
    this.goldenHourEnd = fromJulian(jrise);
    this.goldenHour = fromJulian(jset);
  }
}

class MoonTiming {
  final DateTime date;
  final double lat;
  final double lng;
  bool alwaysUp = false;
  bool alwaysDown = false;
  DateTime? rise;
  DateTime? set;

  MoonTiming(this.date, this.lat, this.lng) {
    const hc = 0.133 * RAD;
    double h0 =
        SunCalc.getMoonPosition(this.date, this.lat, this.lng).altitude - hc;
    var h1 = 0.0;
    var h2 = 0.0;
    var rise = 0.0;
    var set = 0.0;
    var a = 0.0;
    var b = 0.0;
    var xe = 0.0;
    var ye = 0.0;
    var d = 0.0;
    var roots = 0.0;
    var x1 = 0.0;
    var x2 = 0.0;
    var dx = 0.0;

    // go in 2-hour chunks, each time seeing if a 3-point quadratic curve crosses zero (which means rise or set)
    for (var i = 1; i <= 24; i += 2) {
      h1 = SunCalc.getMoonPosition(
                  hoursLater(this.date, i.toDouble()), this.lat, this.lng)
              .altitude -
          hc;
      h2 = SunCalc.getMoonPosition(
                  hoursLater(this.date, i + 1), this.lat, this.lng)
              .altitude -
          hc;

      a = (h0 + h2) / 2 - h1;
      b = (h2 - h0) / 2;
      xe = -b / (2 * a);
      ye = (a * xe + b) * xe + h1;
      d = b * b - 4 * a * h1;
      roots = 0;

      if (d >= 0) {
        dx = math.sqrt(d) / (a.abs() * 2);
        x1 = xe - dx;
        x2 = xe + dx;
        if (x1.abs() <= 1) roots++;
        if (x2.abs() <= 1) roots++;
        if (x1 < -1) x1 = x2;
      }

      if (roots == 1) {
        if (h0 < 0)
          rise = i + x1;
        else
          set = i + x1;
      } else if (roots == 2) {
        rise = i + (ye < 0 ? x2 : x1);
        set = i + (ye < 0 ? x1 : x2);
      }

      if ((rise != 0) && (set != 0)) {
        break;
      }

      h0 = h2;
    }

    if (rise != 0) {
      this.rise = hoursLater(this.date, rise);
    }
    if (set != 0) {
      this.set = hoursLater(this.date, set);
    }

    if ((rise == 0) && (set == 0)) {
      if (ye > 0) {
        this.alwaysUp = true;
      } else {
        this.alwaysDown = true;
      }
    }
  }
}

class CelestialPosition {
  final double azimuth;
  final double altitude;
  final double? distance;
  final double? parallacticAngle;

  CelestialPosition(
      this.azimuth, this.altitude, this.distance, this.parallacticAngle);
}

class MoonIllumination {
  double fraction = 0.0;
  double phase = 0.0;
  double angle = 0.0;

  MoonIllumination(DateTime date) {
    var d = toDays(date);
    var s = sunCoords(d);
    var m = moonCoords(d);

    var sdist = 149598000; // distance from Earth to Sun in km

    var phi = math.acos(math.sin(s["dec"]!) * math.sin(m["dec"]!) +
        math.cos(s["dec"]!) *
            math.cos(m["dec"]! /*!*/) *
            math.cos(s["ra"]! - m["ra"]!));
    var inc =
        math.atan2(sdist * math.sin(phi), m["dist"]! - sdist * math.cos(phi));
    var angle = math.atan2(
        math.cos(s["dec"]!) * math.sin(s["ra"]! - m["ra"]!),
        math.sin(s["dec"]!) * math.cos(m["dec"]!) -
            math.cos(s["dec"]!) *
                math.sin(m["dec"]!) *
                math.cos(s["ra"]! - m["ra"]!));

    this.fraction = (1 + math.cos(inc)) / 2;
    this.phase = 0.5 + 0.5 * inc * (angle < 0 ? -1 : 1) / PI;
    this.angle = angle;
  }
}

double toJulian(DateTime date) {
  return date.difference(julianEpoch).inSeconds / Duration.secondsPerDay;
}

DateTime? fromJulian(double j) {
  if (j.isNaN)
    return null;
  else
    return julianEpoch
        .add(Duration(milliseconds: (j * Duration.millisecondsPerDay).floor()));
}

double toDays(DateTime date) {
  return toJulian(date) - J2000;
}

// general calculations for position
double rightAscension(double l, double b) {
  return math.atan2(
      math.sin(l) * math.cos(E) - math.tan(b) * math.sin(E), math.cos(l));
}

double declination(double l, double b) {
  return math.asin(
      math.sin(b) * math.cos(E) + math.cos(b) * math.sin(E) * math.sin(l));
}

double azimuth(double H, double phi, double dec) {
  return math.atan2(
      math.sin(H), math.cos(H) * math.sin(phi) - math.tan(dec) * math.cos(phi));
}

double altitude(double H, double phi, double dec) {
  return math.asin(math.sin(phi) * math.sin(dec) +
      math.cos(phi) * math.cos(dec) * math.cos(H));
}

double siderealTime(double d, double lw) {
  return RAD * (280.16 + 360.9856235 * d) - lw;
}

double astroRefraction(double h) {
  if (h < 0) {
    // the following formula works for positive altitudes only.
    h = 0; // if h = -0.08901179 a div/0 would occur.
  }
  // formula 16.4 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
  // 1.02 / tan(h + 10.26 / (h + 5.10)) h in degrees, result in arc minutes -> converted to rad:
  return 0.0002967 / math.tan(h + 0.00312536 / (h + 0.08901179));
}

// general sun calculations
double solarMeanAnomaly(double d) {
  return RAD * (357.5291 + 0.98560028 * d);
}

double equationOfCenter(double M) {
  var firstFactor = 1.9148 * math.sin(M);
  var secondFactor = 0.02 * math.sin(2 * M);
  var thirdFactor = 0.0003 * math.sin(3 * M);

  return RAD * (firstFactor + secondFactor + thirdFactor);
}

double eclipticLongitude(double M) {
  var C = equationOfCenter(M);
  var P = RAD * 102.9372; // perihelion of the Earth

  return M + C + P + PI;
}

Map<String, double> sunCoords(double d) {
  var M = solarMeanAnomaly(d);
  var L = eclipticLongitude(M);

  return {"dec": declination(L, 0), "ra": rightAscension(L, 0)};
}

double julianCycle(double d, double lw) {
  return (d - J0 - lw / (2 * PI)).round().toDouble();
}

double approxTransit(double ht, double lw, double n) {
  return J0 + (ht + lw) / (2 * PI) + n;
}

double solarTransitJ(ds, M, L) {
  return J2000 + ds + 0.0053 * math.sin(M) - 0.0069 * math.sin(2 * L);
}

double hourAngle(h, phi, d) {
  return math.acos((math.sin(h) - math.sin(phi) * math.sin(d)) /
      (math.cos(phi) * math.cos(d)));
}

double getSetJ(h, lw, phi, dec, n, M, L) {
  var w = hourAngle(h, phi, dec);
  var a = approxTransit(w, lw, n);

  return solarTransitJ(a, M, L);
}

DateTime hoursLater(DateTime date, double h) {
  var ms = h * 60 * 60 * 1000;
  return date.add(new Duration(milliseconds: ms.toInt()));
}

// moon calculations, based on http://aa.quae.nl/en/reken/hemelpositie.html formulas
Map<String, double> moonCoords(double d) {
  var L = RAD * (218.316 + 13.176396 * d);
  var M = RAD * (134.963 + 13.064993 * d);
  var F = RAD * (93.272 + 13.229350 * d);

  var l = L + RAD * 6.289 * math.sin(M);
  var b = RAD * 5.128 * math.sin(F);
  var dt = 385001 - 20905 * math.cos(M);

  return {"ra": rightAscension(l, b), "dec": declination(l, b), "dist": dt};
}

class SunCalc {
  // calculates sun position for a given date and latitude/longitude
  static CelestialPosition getPosition(DateTime date, double lat, double lng) {
    var lw = RAD * -lng;
    var phi = RAD * lat;
    var d = toDays(date);

    var c = sunCoords(d);
    var H = siderealTime(d, lw) - c["ra"]!;

    return CelestialPosition(
        azimuth(H, phi, c["dec"]!), altitude(H, phi, c["dec"]!), null, null);
  }

  static CelestialPosition getSunPosition(
      DateTime date, double lat, double lng) {
    return SunCalc.getPosition(date, lat, lng);
  }

  static SunTiming getTimes(DateTime date, double lat, double lng) {
    return SunTiming(date, lat, lng);
  }

  static CelestialPosition getMoonPosition(
      DateTime date, double lat, double lng) {
    var lw = RAD * -lng;
    var phi = RAD * lat;
    var d = toDays(date);

    var c = moonCoords(d);
    var H = siderealTime(d, lw) - c["ra"]!;
    var h = altitude(H, phi, c["dec"]!);
    // formula 14.1 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
    var pa = math.atan2(
        math.sin(H),
        math.tan(phi) * math.cos(c["dec"]!) -
            math.sin(c["dec"]!) * math.cos(H));

    h = h + astroRefraction(h); // altitude correction for refraction

    return CelestialPosition(azimuth(H, phi, c["dec"]!), h, c["dist"]!, pa);
  }

  static MoonIllumination getMoonIllumination(DateTime date) {
    return MoonIllumination(date);
  }

  static MoonTiming getMoonTimes(DateTime date, double lat, double lng,
      [bool inUtc = true]) {
    final calcUsingDate = inUtc
        ? new DateTime.utc(date.year, date.month, date.day, 0, 0, 0)
        : new DateTime(date.year, date.month, date.day, 0, 0, 0);
    return MoonTiming(calcUsingDate, lat, lng);
  }
}
