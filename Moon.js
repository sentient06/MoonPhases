//  Created by Giancarlo Mariot in January 2014.
//  Copyright (c) 2012 Giancarlo Mariot. All rights reserved.
//  Licensed under BSD 2-Clause.

/// ----------------
/// Move these away:

function isset(variable) {
    return (typeof(variable) !== 'undefined' && variable !== null);
}

/// ----------------


/// http://stackoverflow.com/questions/2526815/moon-lunar-phase-algorithm
/// http://bazaar.launchpad.net/~keturn/py-moon-phase/trunk/annotate/head:/moon.py

// Based on code by John Walker (http://www.fourmilab.ch/)
// ported to Python by Kevin Turner <acapnotic@twistedmatrix.com>
// on June 6, 2001 (JDN 2452066.52491), under a full moon.

/// Artificial classes:

/**
 * Artificial DateTime class to support proper date.
 * References:
 * http://home.online.no/~pjacklam/matlab/software/util/timeutil/
 * https://en.wikipedia.org/wiki/Julian_date#cite_ref-24
 */
function DateTime(newYear, newMonth, newDay, newHour, newMinute, newSecond, newMicrosecond) {

    "use strict";

    var _this = this;
    var today = new Date();

    // Class' variables
    _this.year   = isset(newYear)   ? newYear   : today.getUTCFullYear();
    _this.month  = isset(newMonth)  ? newMonth  : today.getUTCMonth() + 1;
    _this.date   = isset(newDay)    ? newDay    : today.getUTCDate();
    _this.hour   = isset(newHour)   ? newHour   : today.getUTCHours();
    _this.minute = isset(newMinute) ? newMinute : today.getUTCMinutes();
    _this.second = isset(newSecond) ? newSecond : today.getUTCSeconds();

    _this.DateTimeFromJDN = function(JDN) {

        var y = 4716;
        var j = 1401;
        var m = 2;
        var n = 12;
        var r = 4;
        var p = 1461;
        var v = 3;
        var u = 5;
        var s = 153;
        var w = 2;
        var B = 274277;
        var C = -38;
        var f = JDN + j + (((4 * JDN + B)/146097) * 3)/4 + C;
        var e = r * f + v;
        var g = (e % p)/r;
        var h = u * g + w;
        var D = ((h % s))/u + 1;
        var M = ((h/s + m) % n) + 1;
        var Y = e/p - y + (n + m - M)/n;

        _this.date   = Math.floor(D);
        _this.month  = Math.floor(M);
        _this.year   = Math.floor(Y);
        _this.hour   = 12;
        _this.minute = 0;
        _this.second = 0;

        return _this;
    };

    _this.julianDayNumber = function(newYear, newMonth, newDay, newHour, newMinute, newSecond) {
        var year   = isset(newYear)   ? newYear   : _this.year;
        var month  = isset(newMonth)  ? newMonth  : _this.month;
        var date   = isset(newDay)    ? newDay    : _this.date;
        var hour   = isset(newHour)   ? newHour   : _this.hour;
        var minute = isset(newMinute) ? newMinute : _this.minute;
        var second = isset(newSecond) ? newSecond : _this.second;

        var a = Math.floor((14-month)/12);
        var y = year + 4800 - a;
        var m = month + 12 * a - 3;

        var jdnt =  date + Math.floor((153 * m + 2) / 5) + 365 * y +
                    Math.floor(y/4) - Math.floor(y/100) +
                    Math.floor(y/400) - 32045 + (
                       (hour-12)/24 + minute/1440 + second/86400
                    );
        var jdntInt = Math.floor(jdnt * 1000000);
        return (jdntInt/1000000);

    };

    _this.toString = function(format) {
        return _this.date + "/" + _this.month + "/" + _this.year + " " + _this.hour + ":" + _this.minute + ":" + _this.second;
    };

    _this.RelativeDateTime = function(days) {
        var JDN = _this.julianDayNumber(_this.year, _this.month, _this.date, _this.hour, _this.minute, _this.second);
        JDN += days;
        return _this.DateTimeFromJDN(JDN);
    };

    _this.jdn = _this.julianDayNumber(_this.year, _this.month, _this.date);

}


/*
 * Functions to find the phase of the moon.
 *
 * Ported from "A Moon for the Sun" (aka moontool.c), a program by the
 * venerable John Walker.  He used algoritms from "Practical Astronomy
 * With Your Calculator" by Peter Duffett-Smith, Second Edition.
 *
 * For the full history of the code, as well as references to other
 * reading material and other entertainments, please refer to John
 * Walker's website,
 * http://www.fourmilab.ch/
 * (Look under the Science/Astronomy and Space heading.)
 *
 * The functions of primary interest provided by this module are phase(),
 * which gives you a variety of data on the status of the moon for a
 * given date; and phase_hunt(), which given a date, finds the dates of
 * the nearest full moon, new moon, etc.
 */

// Precision used when describing the moon's phase in textual format,
// in phase_string().

const PRECISION = 0.05;
const NEW       = 0 / 4.0;
const FIRST     = 1 / 4.0;
const FULL      = 2 / 4.0;
const LAST      = 3 / 4.0;
const NEXTNEW   = 4 / 4.0;

/**
 * Little handy mathematical functions.
 */
function Maths() {
    "use strict";
    var _this = this;
    _this.fixangle = function(a) {
        return a - 360.0 * Math.floor(a/360.0);
    };
    _this.torad = function(d) {
        return d * Math.PI / 180.0;
    };
    _this.todeg = function(r) {
        return r * 180.0 / Math.PI;
    };
    _this.dsin = function(d) {
        return Math.sin(_this.torad(d));
    };
    _this.dcos = function(d) {
        return Math.cos(_this.torad(d));
    };
    _this.compareLessThanElement = function (element, index, array) {
        return (this < element);
    };
}

function AstronomicalConstants() {

    "use strict";

    var _this = this;

    // JDN stands for Julian Day Number
    // Angles here are in degrees

    // 1980 January 0.0 in JDN
    // XXX: DateTime(1980).jdn yields 2444239.5 -- which one is right?
    _this.epoch = 2444239.5; //2444238.5;

    // Ecliptic longitude of the Sun at epoch 1980.0
    _this.ecliptic_longitude_epoch = 278.833540;

    // Ecliptic longitude of the Sun at perigee
    _this.ecliptic_longitude_perigee = 282.596403;

    // Eccentricity of Earth's orbit
    _this.eccentricity = 0.016718;

    // Semi-major axis of Earth's orbit, in kilometers
    _this.sun_smaxis = 1.49585e8;

    // Sun's angular size, in degrees, at semi-major axis distance
    _this.sun_angular_size_smaxis = 0.533128;

    // // Elements of the Moon's orbit, epoch 1980.0

    // Moon's mean longitude at the epoch
    _this.moon_mean_longitude_epoch = 64.975464;
    // Mean longitude of the perigee at the epoch
    _this.moon_mean_perigee_epoch = 349.383063;

    // Mean longitude of the node at the epoch
    _this.node_mean_longitude_epoch = 151.950429;

    // Inclination of the Moon's orbit
    _this.moon_inclination = 5.145396;

    // Eccentricity of the Moon's orbit
    _this.moon_eccentricity = 0.054900;

    // Moon's angular size at distance a from Earth
    _this.moon_angular_size = 0.5181;

    // Semi-mojor axis of the Moon's orbit, in kilometers
    _this.moon_smaxis = 384401.0;
    // Parallax at a distance a from Earth
    _this.moon_parallax = 0.9507;

    // Synodic month (new Moon to new Moon), in days
    _this.synodic_month = 29.53058868;

    // Base date for E. W. Brown's numbered series of lunations (1923 January 16)
    _this.lunations_base = 2423436.0;

    // // Properties of the Earth

    _this.earth_radius = 6378.16;

}

/**
 * I describe the phase of the moon.
 * 
 * I have the following properties:
 *     date - a DateTime instance
 *     phase - my phase, in the range 0.0 .. 1.0
 *     phase_text - a string describing my phase
 *     illuminated - the percentage of the face of the moon illuminated
 *     angular_diameter - as seen from Earth, in degrees.
 *     sun_angular_diameter - as seen from Earth, in degrees.
 * 
 *     new_date - the date of the most recent new moon
 *     q1_date - the date the moon reaches 1st quarter in this cycle
 *     full_date - the date of the full moon in this cycle
 *     q3_date - the date the moon reaches 3rd quarter in this cycle
 *     nextnew_date - the date of the next new moon
 */
function MoonPhase(date) {

    "use strict";

    var _self = this;
    _self.c = new AstronomicalConstants();

    /**
     * MoonPhase constructor.
     * Give me a date, as either a Julian Day Number or a DateTime
     * object.
     */
    _self.init = function(date) {
        if (date instanceof DateTime === false) {
            var dateTime = new DateTime();
            _self.date = dateTime.DateTimeFromJDN(date);
        } else {
            _self.date = date;
        }
        var thePhase = _self.phase(_self.date);
        _self.phase_text = _self.phase_string(thePhase.phase);
    };

    /// dead method?
    _self.getattr = function(self, a) {
        var stuff = ['new_date', 'q1_date', 'full_date', 'q3_date', 'nextnew_date'];
        if (stuff.indexOf(a) > -1) {
            //// (self.new_date, self.q1_date, self.full_date, self.q3_date, self.nextnew_date) = phase_hunt(self.date);
            return _self.getattr(self,a);
        }
        console.error("Error! getattr()");
    };

    _self.repr = function(self) {
        if (self.date instanceof int) {
            jdn = self.date;
        } else {
            jdn = self.date.jdn;
        }
        return "<MoonPhase(" + jdn + ")>";
    };

    _self.str = function(self) {
        if (self.date instanceof int) {
            d = DateTime.DateTimeFromJDN(self.date);
        } else {
            d = self.date;
        }
        var s = "MoonPhase for " + d.toString() + ", " + self.phase_text + " (" + (self.illuminated * 100) + " illuminated)";
        return s;
    }

    _self.phase_string = function(p) {
        var phase_keys = [
            NEW     + PRECISION,
            FIRST   - PRECISION,
            FIRST   + PRECISION,
            FULL    - PRECISION,
            FULL    + PRECISION,
            LAST    - PRECISION,
            LAST    + PRECISION,
            NEXTNEW - PRECISION,
            NEXTNEW + PRECISION
        ];
        var phase_strings = [
            "new",
            "waxing crescent (new)",
            "first quarter",
            "waxing gibbous (1st 1/4)",
            "full",
            "waning gibbous (full)",
            "last quarter (3rd 1/4)",
            "waning crescent (3rd 1/4)",
            "new"
        ];
        var maths = new Maths();
        var phaseIndex = phase_keys.findIndex(maths.compareLessThanElement, p);
        var moonPhase = phase_strings[phaseIndex];
        return moonPhase;
    };

    /**
     * Calculate phase of moon as a fraction:
     *
     * The argument is the time for which the phase is requested,
     * expressed in either a DateTime or by Julian Day Number.
     *
     * Returns a dictionary containing the terminator phase angle as a
     * percentage of a full circle (i.e., 0 to 1), the illuminated
     * fraction of the Moon's disc, the Moon's age in days and fraction,
     * the distance of the Moon from the centre of the Earth, and the
     * angular diameter subtended by the Moon as seen by an observer at
     * the centre of the Earth.
     */
    _self.phase = function(phase_date) {
        var c = _self.c;
        var maths = new Maths();

        if (!isset(phase_date)) {
            phase_date = new DateTime();
        }

        // Calculation of the Sun's position
        var day;
        // date within the epoch
        if ('jdn' in phase_date) {
            day = phase_date.jdn - c.epoch;
        } else {
            day = phase_date - c.epoch;
        }

        // Mean anomaly of the Sun
        var N = maths.fixangle((360/365.2422) * day);
        // Convert from perigee coordinates to epoch 1980
        var M = maths.fixangle(N + c.ecliptic_longitude_epoch - c.ecliptic_longitude_perigee);

        // Solve Kepler's equation
        var Ec = null;
        Ec = _self.kepler(M, c.eccentricity);
        Ec = Math.sqrt((1 + c.eccentricity) / (1 - c.eccentricity)) * Math.tan(Ec/2.0);
        // True anomaly
        Ec = 2 * maths.todeg(Math.atan(Ec));
        // Suns's geometric ecliptic longuitude
        var lambda_sun = maths.fixangle(Ec + c.ecliptic_longitude_perigee);

        // Orbital distance factor
        var F = ((1 + c.eccentricity * Math.cos(maths.torad(Ec))) / (1 - Math.pow(c.eccentricity,2)));

        // Distance to Sun in km
        var sun_dist = c.sun_smaxis / F;
        var sun_angular_diameter = F * c.sun_angular_size_smaxis;

        //#######
        //
        // Calculation of the Moon's position

        // Moon's mean longitude
        var moon_longitude = maths.fixangle(13.1763966 * day + c.moon_mean_longitude_epoch);

        // Moon's mean anomaly
        var MM = maths.fixangle(moon_longitude - 0.1114041 * day - c.moon_mean_perigee_epoch);

        // Moon's ascending node mean longitude
        // MN = maths.fixangle(c.node_mean_longitude_epoch - 0.0529539 * day)

        var evection = 1.2739 * Math.sin(maths.torad(2*(moon_longitude - lambda_sun) - MM));

        // Annual equation
        var annual_eq = 0.1858 * Math.sin(maths.torad(M));

        // Correction term
        var A3 = 0.37 * Math.sin(maths.torad(M));

        var MmP = MM + evection - annual_eq - A3;

        // Correction for the equation of the centre
        var mEc = 6.2886 * Math.sin(maths.torad(MmP));

        // Another correction term
        var A4 = 0.214 * Math.sin(maths.torad(2 * MmP));

        // Corrected longitude
        var lP = moon_longitude + evection + mEc - annual_eq + A4;

        // Variation
        var variation = 0.6583 * Math.sin(maths.torad(2*(lP - lambda_sun)));

        // True longitude
        var lPP = lP + variation;

        //
        // Calculation of the Moon's inclination
        // unused for phase calculation.

        // Corrected longitude of the node
        // NP = MN - 0.16 * sin(torad(M))

        // Y inclination coordinate
        // y = sin(torad(lPP - NP)) * cos(torad(c.moon_inclination))

        // X inclination coordinate
        // x = cos(torad(lPP - NP))

        // Ecliptic longitude (unused?)
        // lambda_moon = todeg(atan2(y,x)) + NP

        // Ecliptic latitude (unused?)
        // BetaM = todeg(asin(sin(torad(lPP - NP)) * sin(torad(c.moon_inclination))))

        //######
        //
        // Calculation of the phase of the Moon

        // Age of the Moon, in degrees
        var moon_age = lPP - lambda_sun;

        // Phase of the Moon
        var moon_phase = (1 - Math.cos(maths.torad(moon_age))) / 2.0;

        // Calculate distance of Moon from the centre of the Earth
        var moon_dist = (c.moon_smaxis * (1 - Math.pow(c.moon_eccentricity,2))) /
                        (1 + c.moon_eccentricity * Math.cos(maths.torad(MmP + mEc)));

        // Calculate Moon's angular diameter
        var moon_diam_frac = moon_dist / c.moon_smaxis;
        var moon_angular_diameter = c.moon_angular_size / moon_diam_frac;

        // Calculate Moon's parallax (unused?)
        // moon_parallax = c.moon_parallax / moon_diam_frac

        var res = {
                           'phase': maths.fixangle(moon_age) / 360.0,
                     'illuminated': moon_phase,
                             'age': c.synodic_month * maths.fixangle(moon_age) / 360.0 ,
                        'distance': moon_dist,
                'angular_diameter': moon_angular_diameter,
                    'sun_distance': sun_dist,
            'sun_angular_diameter': sun_angular_diameter
        };

        return res;
    };

    /**
     * Find time of phases of the moon which surround the current date.
     *
     * Five phases are found, starting and ending with the new moons
     * which bound the current lunation.
     */
    _self.phase_hunt = function(sdate) {

        var c = _self.c;
        var dateTime = new DateTime();

        if (!isset(sdate)) {
            sdate = new DateTime();
        }

        var adate = sdate.RelativeDateTime(-45); /// Test!
        var k1 = Math.floor((adate.year + ((adate.month - 1) * (1.0/12.0)) - 1900) * 12.3685);
        var nt1 = _self.meanphase(adate.julianDayNumber(), k1);
        adate = nt1;
        sdate = sdate.jdn;

        var k2, nt2;

        while(1){
            adate = adate + c.synodic_month;
            k2 = k1 + 1;
            nt2 = _self.meanphase(adate,k2);
            if (nt1 <= sdate < nt2) {
                break;
            }
            nt1 = nt2;
            k1 = k2;
        }

        var phases = [
            _self.truephase(k1, 0/4.0),
            _self.truephase(k1, 1/4.0),
            _self.truephase(k1, 2/4.0),
            _self.truephase(k1, 3/4.0),
            _self.truephase(k2, 0/4.0)
        ];

        return phases;
    };


    /**
     * Calculate mean new moon
     *
     * Calculates  time  of  the mean new Moon for a given
     * base date.  This argument K to this function is the
     * precomputed synodic month index, given by:
     *     K = (year - 1900) * 12.3685
     * where year is expressed as a year and fractional year.
     *
     * @param number sdate - it is a Julian date
     * @param number k
     * @return number
     */
    _self.meanphase = function (sdate, k) {
        // Time in Julian centuries from 1900 January 0.5

        // JD 2415020.000000 is
        // CE 1899 December 31 12:00:00.0 UT Sunday
        var c = _self.c;
        var maths = new Maths();

        var time1 = ( sdate - 2415020.0 ) / 36525;
        var time2 = time1 * time1;
        var time3 = time2 * time1;

        return 2415020.75933 + c.synodic_month * k
            + 0.0001178 * time2
            - 0.000000155 * time3
            + 0.00033 * maths.dsin(166.56 + 132.87 * time1 - 0.009173 * time2);
    };

    /**
     * Given a K value used to determine the mean phase of the new
     * moon, and a phase selector (0.0, 0.25, 0.5, 0.75), obtain the
     * true, corrected phase time.
     */
    _self.truephase = function(k, tphase) {

        var c = _self.c;
        var maths = new Maths();

        var apcor = false;

        // add phase to new moon time
        k = k + tphase;
        // Time in Julian centuries from 1900 January 0.5
        var t = k / 1236.85;

        var t2 = t * t;
        var t3 = t2 * t;

        // Mean time of phase
        var pt = 2415020.75933 + c.synodic_month * k + 0.0001178 * t2 -
            0.000000155 * t3 + 0.00033 * maths.dsin(166.56 + 132.87 * t - 0.009173 * t2);

        // Sun's mean anomaly
        var m = 359.2242 + 29.10535608 * k - 0.0000333 * t2 - 0.00000347 * t3;

        // Moon's mean anomaly
        var mprime = 306.0253 + 385.81691806 * k + 0.0107306 * t2 + 0.00001236 * t3;

        // Moon's argument of latitude
        var f = 21.2964 + 390.67050646 * k - 0.0016528 * t2 - 0.00000239 * t3;

        if ((tphase < 0.01) || (Math.abs(tphase - 0.5) < 0.01)) {

            // Corrections for New and Full Moon
            pt = pt + (
                (0.1734 - 0.000393 * t) * maths.dsin(m)
                + 0.0021 * maths.dsin(2 * m)
                - 0.4068 * maths.dsin(mprime)
                + 0.0161 * maths.dsin(2 * mprime)
                - 0.0004 * maths.dsin(3 * mprime)
                + 0.0104 * maths.dsin(2 * f)
                - 0.0051 * maths.dsin(m + mprime)
                - 0.0074 * maths.dsin(m - mprime)
                + 0.0004 * maths.dsin(2 * f + m)
                - 0.0004 * maths.dsin(2 * f - m)
                - 0.0006 * maths.dsin(2 * f + mprime)
                + 0.0010 * maths.dsin(2 * f - mprime)
                + 0.0005 * maths.dsin(m + 2 * mprime)
            );

            apcor = true;

        } else if ((Math.abs(tphase - 0.25) < 0.01) || (Math.abs(tphase - 0.75) < 0.01)) {

            pt = pt + (
                (0.1721 - 0.0004 * t) * maths.dsin(m)
                + 0.0021 * maths.dsin(2 * m)
                - 0.6280 * maths.dsin(mprime)
                + 0.0089 * maths.dsin(2 * mprime)
                - 0.0004 * maths.dsin(3 * mprime)
                + 0.0079 * maths.dsin(2 * f)
                - 0.0119 * maths.dsin(m + mprime)
                - 0.0047 * maths.dsin(m - mprime)
                + 0.0003 * maths.dsin(2 * f + m)
                - 0.0004 * maths.dsin(2 * f - m)
                - 0.0006 * maths.dsin(2 * f + mprime)
                + 0.0021 * maths.dsin(2 * f - mprime)
                + 0.0003 * maths.dsin(m + 2 * mprime)
                + 0.0004 * maths.dsin(m - 2 * mprime)
                - 0.0003 * maths.dsin(2 * m + mprime)
            );

            if (tphase < 0.5) {
                // First quarter correction
                pt = pt + 0.0028 - 0.0004 * maths.dcos(m) + 0.0003 * maths.dcos(mprime);
            } else {
                //  Last quarter correction
                pt = pt + -0.0028 + 0.0004 * maths.dcos(m) - 0.0003 * maths.dcos(mprime);
            }

            apcor = true;
        }

        if (!apcor) {
            console.error("TRUEPHASE called with invalid phase selector", tphase);
        }

        var dateTime = new DateTime();

        return dateTime.DateTimeFromJDN(pt);
    }

    /**
     * Solve the equation of Kepler.
     */
    _self.kepler = function(m, ecc) {

        var maths = new Maths();
        var epsilon = Math.pow(1, -6);
        var e = m = maths.torad(m);
        var delta = null;
        do {
            delta = e - ecc * Math.sin(e) - m;
            e -= delta / ( 1 - ecc * Math.cos(e) );
        } while (Math.abs(delta) > epsilon);
        return e;
    }

    _self.init(date);

}

/// Tests

(function(){
    "use strict";
    var doc = document;
    doc.addEventListener('DOMContentLoaded', function() {


    var today = new Date();
    var aDate = new DateTime(today.getUTCFullYear(), today.getUTCMonth() + 1, today.getUTCDate(), 0, 0, 0);

    var moon = new MoonPhase();
    var hunted = moon.phase_hunt(aDate);
    for (var i = 0; i < hunted.length; i++) {
        console.log(hunted[i].toString());
    };
    

    // console.log(moon);
    // console.log(aDate.toString());
    console.log(moon.phase_text);

    },false);
})();
