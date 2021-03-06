#ifndef __ODAS_SYSTEM_ACORR2PITCH
#define __ODAS_SYSTEM_ACORR2PITCH

   /**
    * \file     acorr2pitch.h
    * \author   François Grondin <francois.grondin2@usherbrooke.ca>
    * \version  2.0
    * \date     2018-03-18
    * \copyright
    *
    * This program is free software: you can redistribute it and/or modify
    * it under the terms of the GNU General Public License as published by
    * the Free Software Foundation, either version 3 of the License, or
    * (at your option) any later version.
    *
    * This program is distributed in the hope that it will be useful,
    * but WITHOUT ANY WARRANTY; without even the implied warranty of
    * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    * GNU General Public License for more details.
    * 
    * You should have received a copy of the GNU General Public License
    * along with this program.  If not, see <http://www.gnu.org/licenses/>.
    *
    */

    #include <signal/acorr.h>
    #include <signal/pitch.h>

    #include <math.h>
    // sd included for improved classification
    #include <signal/track.h>


    typedef struct acorr2pitch_obj {

        unsigned int nSignals;
        unsigned int halfFrameSize;
        unsigned int frameSize;
        unsigned int winSize;

    } acorr2pitch_obj;

    acorr2pitch_obj * acorr2pitch_construct_zero(const unsigned int nSignals, const unsigned int halfFrameSize, const unsigned int winSize);

    void acorr2pitch_destroy(acorr2pitch_obj * obj);

    void acorr2pitch_process(acorr2pitch_obj * obj, const acorrs_obj * acorrs, const tracks_obj * tracks, pitches_obj * pitches);

#endif
