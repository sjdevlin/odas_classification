
/**
    * \file     pitch2category.c
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

#include <system/pitch2category.h>

pitch2category_obj *pitch2category_construct_zero(const unsigned int nSeps, const float tauMin, const float tauMax, const float deltaTauMax, const float alpha, const float gamma, const float phiMin, const float r0, const float activityMin, const float acorrMin, const unsigned int classificationPeriod)
{

    pitch2category_obj *obj;

    unsigned int iSep;

    obj = (pitch2category_obj *)malloc(sizeof(pitch2category_obj));

    obj->nSeps = nSeps;

    obj->tauMin = tauMin;
    obj->tauMax = tauMax;
    obj->deltaTauMax = deltaTauMax;
    obj->alpha = alpha;
    obj->gamma = gamma;
    obj->phiMin = phiMin;
    obj->r0 = r0;
    obj->activityMin = activityMin;
    obj->acorrMin = acorrMin;

    //sd changed to improve classification
    obj->classificationPeriod = classificationPeriod;

    obj->tausNow = (float *)malloc(sizeof(float) * nSeps);
    memset(obj->tausNow, 0x00, sizeof(float) * nSeps);

    obj->tausPrev = (float *)malloc(sizeof(float) * nSeps);
    memset(obj->tausPrev, 0x00, sizeof(float) * nSeps);

    obj->deltaTausNow = (float *)malloc(sizeof(float) * nSeps);
    memset(obj->deltaTausNow, 0x00, sizeof(float) * nSeps);

    obj->deltaTausPrev = (float *)malloc(sizeof(float) * nSeps);
    memset(obj->deltaTausPrev, 0x00, sizeof(float) * nSeps);

    obj->phisNow = (float *)malloc(sizeof(float) * nSeps);
    memset(obj->phisNow, 0x00, sizeof(float) * nSeps);

    obj->phisPrev = (float *)malloc(sizeof(float) * nSeps);
    memset(obj->phisPrev, 0x00, sizeof(float) * nSeps);

    obj->vs = (float *)malloc(sizeof(float) * nSeps);
    memset(obj->vs, 0x00, sizeof(float) * nSeps);

    obj->rs = (float *)malloc(sizeof(float) * nSeps);
    memset(obj->rs, 0x00, sizeof(float) * nSeps);

    obj->categories = (char *)malloc(sizeof(char) * nSeps);
    memset(obj->categories, 0x00, sizeof(char) * nSeps);

    // sd changed to improve classification
    obj->processingTime = (unsigned int *)malloc(sizeof(unsigned int) * nSeps);
    memset(obj->processingTime, 0x00, sizeof(unsigned int) * nSeps);

    obj->numPitchValues = (unsigned int *)malloc(sizeof(unsigned int) * nSeps);
    memset(obj->numPitchValues, 0x00, sizeof(unsigned int) * nSeps);

    obj->pitchArray = (unsigned int **)malloc(sizeof(unsigned int *) * nSeps);
    obj->activityArray = (unsigned int **)malloc(sizeof(unsigned int *) * nSeps);

    for (iSep = 0; iSep < nSeps; iSep++)
    {
        obj->pitchArray[iSep] = (unsigned int *)malloc(sizeof(unsigned int) * classificationPeriod);
        memset(obj->pitchArray[iSep], 0x00, sizeof(unsigned int) * classificationPeriod);
        obj->activityArray[iSep] = (unsigned int *)malloc(sizeof(unsigned int) * classificationPeriod);
        memset(obj->activityArray[iSep], 0x00, sizeof(unsigned int) * classificationPeriod);
    }

    obj->pitchTotal = (unsigned int *)malloc(sizeof(unsigned int) * nSeps);
    memset(obj->pitchTotal, 0x00, sizeof(unsigned int) * nSeps);
    obj->activityTotal = (unsigned int *)malloc(sizeof(unsigned int) * nSeps);
    memset(obj->activityTotal, 0x00, sizeof(unsigned int) * nSeps);
    obj->harmonicAcorrTotal = (float *)malloc(sizeof(float) * nSeps);
    memset(obj->harmonicAcorrTotal, 0x00, sizeof(float) * nSeps);

    return obj;
}

void pitch2category_destroy(pitch2category_obj *obj)
{

    free((void *)obj->tausNow);
    free((void *)obj->tausPrev);
    free((void *)obj->deltaTausNow);
    free((void *)obj->deltaTausPrev);
    free((void *)obj->phisNow);
    free((void *)obj->phisPrev);
    free((void *)obj->vs);
    free((void *)obj->rs);
    free((void *)obj->categories);
    // sd changed to improve classification
    free((void *)obj->pitchArray);
    free((void *)obj->pitchTotal);
    free((void *)obj->harmonicAcorrTotal);
    free((void *)obj->numPitchValues);
    free((void *)obj->activityArray);
    free((void *)obj->activityTotal);
    free((void *)obj->processingTime);

    free((void *)obj);
}

void pitch2category_process(pitch2category_obj *obj, const pitches_obj *pitches, const tracks_obj *tracks, categories_obj *categories)
{

    unsigned int iSep;
    float deltaPitch;

    int i;
    float amplitudeMean, pitchMean, harmonicAcorrMean;
    float amplitudeDiff, pitchDiff;
    float totalAmpDiffSquared, totalPitchDiffSquared;
    float relAmpVariance, relPitchVariance;

    for (iSep = 0; iSep < obj->nSeps; iSep++)
    {

        if (tracks->ids[iSep] != 0)
        {

            if (obj->processingTime[iSep] < obj->classificationPeriod) // increment so we know how long we have been looking

            {
                // add activity to array
                obj->activityArray[iSep][obj->processingTime[iSep]] = pitches->realRMS[iSep];
                // also calculate activity total for variance analysis
                obj->activityTotal[iSep] += obj->activityArray[iSep][obj->processingTime[iSep]]; // to avoid iterating a second time

                if (pitches->array[iSep] > 0) // skip analysis when no peak found
                {
                    obj->pitchArray[iSep][obj->numPitchValues[iSep]] = pitches->array[iSep];
                    obj->pitchTotal[iSep] += obj->pitchArray[iSep][obj->numPitchValues[iSep]];
                    obj->harmonicAcorrTotal[iSep] += pitches->harmonicAcorr[iSep];
                    ++obj->numPitchValues[iSep];
                }
            }
            else if (obj->processingTime[iSep] == obj->classificationPeriod) // we have reached time to classify as speech
            {

                // calculate mean of activity
                amplitudeMean = (float)obj->activityTotal[iSep] / (float)obj->classificationPeriod;

                if (obj->numPitchValues[iSep] != 0)
                {
                    pitchMean = (float)obj->pitchTotal[iSep] / (float)obj->numPitchValues[iSep]; // consider changing to int for spped 
                    harmonicAcorrMean = obj->harmonicAcorrTotal[iSep] / (float)obj->numPitchValues[iSep];
                }
                else
                {
                    pitchMean = 0.0f;
                    harmonicAcorrMean = 0.0f;
                }

                totalAmpDiffSquared = 0.0f;
                totalPitchDiffSquared = 0.0f;

                // relative variance of amplitide

                for (i = 0; i < obj->classificationPeriod; i++)
                {
//                    printf (" %d,", obj->activityArray[iSep][i]);
                    amplitudeDiff = (float)obj->activityArray[iSep][i] - amplitudeMean;
                    totalAmpDiffSquared += (amplitudeDiff * amplitudeDiff);  // consider changing to powf
                }

                // relative variance of pitch

                for (i = 0; i < obj->numPitchValues[iSep]; i++)
                {
                    pitchDiff = (float)obj->pitchArray[iSep][i] - pitchMean;  //consider changing to int
                    totalPitchDiffSquared += (pitchDiff * pitchDiff);  // consider powf?
                }

                // rel variance is variance / mean ^2

                relAmpVariance = (totalAmpDiffSquared / (float)obj->classificationPeriod) / (amplitudeMean * amplitudeMean);

                if (obj->numPitchValues[iSep] != 0)
                {
                    relPitchVariance = (totalPitchDiffSquared / (float)obj->numPitchValues[iSep]) / (pitchMean * pitchMean);
                }
                else
                {
                    relPitchVariance = 0.0f;
                }

                if (amplitudeMean > 0.3) printf("%llu, %3.2f, %2.2f, %2.2f, %2.2f,%d \n", tracks->ids[iSep], pitchMean, relPitchVariance, relAmpVariance, harmonicAcorrMean, obj->numPitchValues[iSep]);
//                printf("\nT, %llu, tot diff squared, %5.5f, Amp Mean, %5.5f, Amp Total, %d", tracks->ids[iSep], totalAmpDiffSquared ,amplitudeMean, obj->activityTotal[iSep]); 

            }
            ++obj->processingTime[iSep];
        }
        else
        {

            obj->tausNow[iSep] = 0.0f;
            obj->deltaTausNow[iSep] = 0.0f;
            obj->phisNow[iSep] = 0.0f;
            obj->vs[iSep] = 0.0f;
            obj->rs[iSep] = 0.0f;
            obj->categories[iSep] = 0x00;
            // added for amplitude variance analysis
            obj->processingTime[iSep] = 0;
            obj->activityTotal[iSep] = 0;
            obj->numPitchValues[iSep] = 0;
            obj->pitchTotal[iSep] = 0;
            obj->harmonicAcorrTotal[iSep] = 0.0f;
        }
    }
}
