#ifndef __ODAS_SYSTEM_HOP2FRAME
#define __ODAS_SYSTEM_HOP2FRAME

   /**
    * \file     hop2frame.h
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

    #include <stdlib.h>
    #include <string.h>

    #include <signal/hop.h>
    #include <signal/frame.h>


    typedef struct hop2frame_obj {

        unsigned int hopSize;
        unsigned int frameSize;
        unsigned int nSignals;

        float ** array;

    } hop2frame_obj;

    hop2frame_obj * hop2frame_construct_zero(const unsigned int hopSize, const unsigned int frameSize, const unsigned int nSignals);

    void hop2frame_destroy(hop2frame_obj * obj);

    void hop2frame_process(hop2frame_obj * obj, const hops_obj * hops,frames_obj * frames);

#endif
