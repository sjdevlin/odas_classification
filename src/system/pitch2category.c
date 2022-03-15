
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
    obj->rmsArray = (unsigned int **)malloc(sizeof(unsigned int *) * nSeps);

    for (iSep = 0; iSep < nSeps; iSep++)
    {
        obj->pitchArray[iSep] = (unsigned int *)malloc(sizeof(unsigned int) * classificationPeriod);
        memset(obj->pitchArray[iSep], 0x00, sizeof(unsigned int) * classificationPeriod);
        obj->activityArray[iSep] = (unsigned int *)malloc(sizeof(unsigned int) * classificationPeriod);
        memset(obj->activityArray[iSep], 0x00, sizeof(unsigned int) * classificationPeriod);
        obj->rmsArray[iSep] = (unsigned int *)malloc(sizeof(unsigned int) * classificationPeriod);
        memset(obj->rmsArray[iSep], 0x00, sizeof(unsigned int) * classificationPeriod);
    }

    obj->pitchTotal = (unsigned int *)malloc(sizeof(unsigned int) * nSeps);
    memset(obj->pitchTotal, 0x00, sizeof(unsigned int) * nSeps);
    obj->activityTotal = (unsigned int *)malloc(sizeof(unsigned int) * nSeps);
    memset(obj->activityTotal, 0x00, sizeof(unsigned int) * nSeps);
    obj->rmsTotal = (unsigned int *)malloc(sizeof(unsigned int) * nSeps);
    memset(obj->rmsTotal, 0x00, sizeof(unsigned int) * nSeps);
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
    free((void *)obj->rmsArray);
    free((void *)obj->rmsTotal);
    free((void *)obj->processingTime);

    free((void *)obj);
}

void pitch2category_process(pitch2category_obj *obj, const pitches_obj *pitches, const tracks_obj *tracks, categories_obj *categories, const unsigned int iSep)
{
    int i;
    float deltaPitch;
    float rmsMean, activityMean, pitchMean, harmonicAcorrMean;
    float activityDiff, rmsDiff, pitchDiff;
    float totalRmsDiffSquared, totalActivityDiffSquared, totalPitchDiffSquared;
    float relActivityVariance, relRmsVariance, relPitchVariance;

    // add activity to array
    obj->activityArray[iSep][obj->processingTime[iSep]] = tracks->activity[iSep] * 100;
//    printf ("%d ", obj->activityArray[iSep][obj->processingTime[iSep]]);
    // also calculate activity total for variance analysis
    obj->activityTotal[iSep] += obj->activityArray[iSep][obj->processingTime[iSep]]; // to avoid iterating a second time

    // add rms to array
    obj->rmsArray[iSep][obj->processingTime[iSep]] = pitches->realRMS[iSep];
    // also calculate rms total for variance analysis
    obj->rmsTotal[iSep] += obj->rmsArray[iSep][obj->processingTime[iSep]]; // to avoid iterating a second time

    if (pitches->array[iSep] > 0) // skip analysis when no peak found
    {
        obj->pitchArray[iSep][obj->numPitchValues[iSep]] = pitches->array[iSep];
        obj->pitchTotal[iSep] += obj->pitchArray[iSep][obj->numPitchValues[iSep]];
        obj->harmonicAcorrTotal[iSep] += pitches->harmonicAcorr[iSep];
        ++obj->numPitchValues[iSep];
    }

    // when classification period is ended do all the processing

    if (obj->processingTime[iSep] == obj->classificationPeriod) // we have reached time to classify as speech
    {

        // calculate mean of activity
        activityMean = (float)obj->activityTotal[iSep] / (float)obj->classificationPeriod;
        // calculate mean of rms
        rmsMean = (float)obj->rmsTotal[iSep] / (float)obj->classificationPeriod;

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

        totalActivityDiffSquared = 0.0f;
        totalRmsDiffSquared = 0.0f;
        totalPitchDiffSquared = 0.0f;

        // relative variance of amplitide and E (for comparison)

        for (i = 0; i < obj->classificationPeriod; i++)
        {
            //                    printf (" %d,", obj->activityArray[iSep][i]);
            activityDiff = (float)obj->activityArray[iSep][i] - activityMean;
            totalActivityDiffSquared += (activityDiff * activityDiff); // consider changing to powf

            rmsDiff = (float)obj->rmsArray[iSep][i] - rmsMean;
            totalRmsDiffSquared += (rmsDiff * rmsDiff); // consider changing to powf
        }

        // relative variance of pitch

        for (i = 0; i < obj->numPitchValues[iSep]; i++)
        {
            pitchDiff = (float)obj->pitchArray[iSep][i] - pitchMean; //consider changing to int
            totalPitchDiffSquared += (pitchDiff * pitchDiff);        // consider powf?
        }

        // rel variance is variance / mean ^2

       relActivityVariance = (totalActivityDiffSquared / (float)obj->classificationPeriod) / (activityMean * activityMean);
//    printf("Act mean: %2.3f, Act Diff Sq: %2.3f \n" ,  activityMean, totalActivityDiffSquared);
       relRmsVariance = (totalRmsDiffSquared / (float)obj->classificationPeriod) / (rmsMean * rmsMean);
 
        if (obj->numPitchValues[iSep] != 0)
        {
            relPitchVariance = (totalPitchDiffSquared / (float)obj->numPitchValues[iSep]) / (pitchMean * pitchMean);
        }
        else
        {
            relPitchVariance = 0.0f;
        }

        printf("%llu, %3.2f, %2.2f, %2.2f, %2.2f, %2.2f, %d \n", tracks->ids[iSep], pitchMean, relPitchVariance, relRmsVariance, relActivityVariance, harmonicAcorrMean, obj->numPitchValues[iSep]);
        //                printf("\nT, %llu, tot diff squared, %5.5f, Amp Mean, %5.5f, Amp Total, %d", tracks->ids[iSep], totalAmpDiffSquared ,amplitudeMean, obj->activityTotal[iSep]);

        // reset everything to zero
    }

}