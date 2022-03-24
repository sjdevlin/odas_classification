
/**
    * \file     mod_classify.c
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

#include <module/mod_classify.h>

mod_classify_obj *mod_classify_construct(const mod_classify_cfg *mod_classify_config, const msg_hops_cfg *msg_hops_config, const msg_tracks_cfg *msg_tracks_config, const msg_categories_cfg *msg_categories_config)
{

    mod_classify_obj *obj;

    obj = (mod_classify_obj *)malloc(sizeof(mod_classify_obj));

    obj->hop2frame = hop2frame_construct_zero(msg_hops_config->hopSize,
                                              mod_classify_config->frameSize,
                                              msg_hops_config->nChannels);

    obj->frames = frames_construct_zero(msg_hops_config->nChannels,
                                        mod_classify_config->frameSize);

    obj->frame2freq = frame2freq_construct_zero(mod_classify_config->frameSize,
                                                (mod_classify_config->frameSize) / 2 + 1);

    obj->freqs = freqs_construct_zero(msg_hops_config->nChannels,
                                      (mod_classify_config->frameSize) / 2 + 1);

    obj->freq2acorr = freq2acorr_construct_zero(msg_hops_config->nChannels,
                                                (mod_classify_config->frameSize) / 2 + 1);

    obj->acorrs = acorrs_construct_zero(msg_hops_config->nChannels,
                                        (mod_classify_config->frameSize) / 2 + 1);

    obj->acorr2pitch = acorr2pitch_construct_zero(msg_hops_config->nChannels,
                                                  (mod_classify_config->frameSize) / 2 + 1,
                                                  mod_classify_config->winSize);

    obj->pitches = pitches_construct_zero(msg_hops_config->nChannels);

    obj->pitch2category = pitch2category_construct_zero(msg_hops_config->nChannels,
                                                        mod_classify_config->tauMin,
                                                        mod_classify_config->tauMax,
                                                        mod_classify_config->deltaTauMax,
                                                        mod_classify_config->alpha,
                                                        mod_classify_config->gamma,
                                                        mod_classify_config->phiMin,
                                                        mod_classify_config->r0,
                                                        mod_classify_config->activityMin,
                                                        mod_classify_config->acorrMin,
                                                        mod_classify_config->classificationPeriod);

    obj->in1 = (msg_hops_obj *)NULL;
    obj->in2 = (msg_tracks_obj *)NULL;
    obj->out = (msg_categories_obj *)NULL;

    obj->enabled = 0;

    return obj;
}

void mod_classify_destroy(mod_classify_obj *obj)
{

    hop2frame_destroy(obj->hop2frame);
    frames_destroy(obj->frames);
    frame2freq_destroy(obj->frame2freq);
    freqs_destroy(obj->freqs);
    freq2acorr_destroy(obj->freq2acorr);
    acorrs_destroy(obj->acorrs);
    acorr2pitch_destroy(obj->acorr2pitch);
    pitches_destroy(obj->pitches);
    pitch2category_destroy(obj->pitch2category);

    free((void *)obj);
}

int mod_classify_process(mod_classify_obj *obj)
{
    unsigned int iSignal;

    float Yreal;
    float Yimag;
    float Ypwr;

    int rtnValue;

    signed int iBin;
    signed int iWin;
    signed int iSample;

    float peakValue;
    float nextValue;
    float maxValue;

    unsigned int maxIndex;
    char isPeak;
    char peakFound;

    unsigned int peakPrev;
    unsigned int peakNow;
    unsigned int peakNext;

    float RPrev;
    float RNow;
    float RNext;

    float totalRealAmplitudeSquared;
    float meanRealAmplitudeSquared;
    float realRMS;

    if (obj->in1->timeStamp != obj->in2->timeStamp)
    {
        printf("Time stamp mismatch.\n");
        exit(EXIT_FAILURE);
    }

    if (msg_hops_isZero(obj->in1) == 0)
    {

        if (obj->enabled == 1)
        {
            //   Here we need to make sure we only process this once !
            for (iSignal = 0; iSignal < obj->in2->tracks->nTracks; iSignal++)
            {

                if (obj->in2->tracks->ids[iSignal] != 0)
                {
                    if (obj->pitch2category->categories[iSignal] == 0x00)
                    {
//                        printf("%d, ", obj->frames->numHops[iSignal]);

                        memmove(&(obj->hop2frame->array[iSignal][0]),
                                &(obj->hop2frame->array[iSignal][obj->hop2frame->hopSize]),
                                sizeof(float) * (obj->hop2frame->frameSize - obj->hop2frame->hopSize));

                        memcpy(&(obj->hop2frame->array[iSignal][(obj->hop2frame->frameSize - obj->hop2frame->hopSize)]),
                               &(obj->in1->hops->array[iSignal][0]),
                               sizeof(float) * obj->hop2frame->hopSize);

                        memcpy(obj->frames->array[iSignal],
                               obj->hop2frame->array[iSignal],
                               sizeof(float) * obj->hop2frame->frameSize);

                        // We only want to do the following processes every time we get to
                        // half the classification frame size (for each signal) - 

                        if (obj->pitch2category->processingTime[iSignal] <= obj->pitch2category->classificationPeriod)
                        {
//                            printf("\n", obj->frames->numHops[iSignal]);
                            ++obj->frames->numHops[iSignal];
                            if (obj->frames->numHops[iSignal] == obj->frame2freq->halfFrameSize / obj->hop2frame->hopSize) //  this needs to be changed to work out the number
                            {

                                obj->frames->numHops[iSignal] = 0;  // reset hop counter
                                totalRealAmplitudeSquared = 0.0f;

                                for (iSample = 0; iSample < obj->frame2freq->frameSize; iSample++)
                                {
                                    // this is where windowing happens - remove the product with win->array to use rectangular window

                                    obj->frame2freq->frame[iSample] = obj->frame2freq->win->array[iSample] * obj->frames->array[iSignal][iSample];
                                    //                                obj->frame2freq->frame[iSample] = obj->frames->array[iSignal][iSample];
                                    totalRealAmplitudeSquared += obj->frame2freq->frameSize * powf(obj->frames->array[iSignal][iSample], 2);
                                }

                                // rms NEEDED FOR classification
                                meanRealAmplitudeSquared = totalRealAmplitudeSquared / (float)obj->frame2freq->frameSize;
                                obj->pitches->realRMS[iSignal] = sqrt(meanRealAmplitudeSquared);

                                if (obj->in2->tracks->activity[iSignal] > obj->pitch2category->activityMin)
                                {

                                    fft_r2c(obj->frame2freq->fft,
                                            obj->frame2freq->frame,
                                            obj->freqs->array[iSignal]);

                                    memset(obj->freq2acorr->arrayIn, 0x00, sizeof(float) * obj->freq2acorr->halfFrameSize * 2);

                                    for (iBin = 0; iBin < obj->freq2acorr->halfFrameSize; iBin++)
                                    {

                                        Yreal = obj->freqs->array[iSignal][iBin * 2 + 0];
                                        Yimag = obj->freqs->array[iSignal][iBin * 2 + 1];
                                        Ypwr = Yreal * Yreal + Yimag * Yimag;

                                        obj->freq2acorr->arrayIn[iBin * 2 + 0] = Ypwr;
                                        obj->freq2acorr->arrayIn[iBin * 2 + 1] = 0.0f;
                                    }

                                    fft_c2r(obj->freq2acorr->fft, obj->freq2acorr->arrayIn, obj->freq2acorr->arrayOut);
                                    memcpy(obj->acorrs->array[iSignal], obj->freq2acorr->arrayOut, sizeof(float) * obj->freq2acorr->halfFrameSize);

                                    maxValue = -INFINITY;
                                    peakFound = 0x00;

                                    for (iBin = obj->acorr2pitch->winSize; iBin < (obj->acorr2pitch->halfFrameSize - obj->acorr2pitch->winSize); iBin++)
                                    {

                                        peakValue = obj->acorrs->array[iSignal][iBin];
                                        isPeak = 0x01;
                                        //                        if (peakValue > 0.05) (printf("accor %2.2f\n", peakValue));

                                        // SD Changes to improve classification
                                        if (peakValue > obj->pitch2category->acorrMin)
                                        {

                                            for (iWin = -1 * ((signed int)(obj->acorr2pitch->winSize)); iWin <= ((signed int)(obj->acorr2pitch->winSize)); iWin++)
                                            {

                                                iSample = iBin + iWin;
                                                nextValue = obj->acorrs->array[iSignal][iSample];

                                                if ((iWin != 0) && (peakValue <= nextValue))
                                                {

                                                    isPeak = 0x00;
                                                    break;
                                                }
                                            }

                                            if (isPeak == 0x01)
                                            {

                                                if (peakValue > maxValue)
                                                {

                                                    maxValue = peakValue;
                                                    maxIndex = iBin;

                                                    peakFound = 0x01;
                                                }
                                            }
                                        }
                                    }

                                    if (peakFound == 0x01)
                                    {
                                        obj->pitches->array[iSignal] = obj->in1->fS / maxIndex; 

                                        if (maxIndex < obj->acorr2pitch->halfFrameSize)
                                        {
                                            obj->pitches->harmonicAcorr[iSignal] = obj->acorrs->array[iSignal][maxIndex * 2] / obj->acorrs->array[iSignal][maxIndex];
                                        }
                                        else
                                        {
                                            obj->pitches->harmonicAcorr[iSignal] = 0.0f;
                                        }
                                    }
                                    else
                                    {
                                        obj->pitches->array[iSignal] = 0;
                                        obj->pitches->harmonicAcorr[iSignal] = 0.0f;
                                    }
                                }
                                else
                                {
                                    obj->pitches->array[iSignal] = 0;
                                    obj->pitches->harmonicAcorr[iSignal] = 0.0f;
                                }
                                if (obj->pitch2category->processingTime[iSignal] <= obj->pitch2category->classificationPeriod)
                                {
                                    pitch2category_process(obj->pitch2category, obj->pitches, obj->in2->tracks, obj->out->categories, iSignal);
                                    ++obj->pitch2category->processingTime[iSignal];
                                }
                            }
                        }
                    }
                }
                else
                {
                    obj->pitch2category->processingTime[iSignal] = 0;
                    obj->frames->numHops[iSignal] = 0;
                    obj->pitch2category->activityTotal[iSignal] = 0;
                    obj->pitch2category->rmsTotal[iSignal] = 0;
                    obj->pitch2category->numPitchValues[iSignal] = 0;
                    obj->pitch2category->pitchTotal[iSignal] = 0;
                    obj->pitch2category->harmonicAcorrTotal[iSignal] = 0.0f;
                    obj->pitch2category->categories[iSignal] = 0x00;
                }
            }
        }
        else
        {
            categories_zero(obj->out->categories);
        }

        obj->out->timeStamp = obj->in1->timeStamp;
        rtnValue = 0;
    }
    else
    {

        msg_categories_zero(obj->out);
        rtnValue = -1;
    }

    return rtnValue;
}

void mod_classify_connect(mod_classify_obj *obj, msg_hops_obj *in1, msg_tracks_obj *in2, msg_categories_obj *out)
{

    obj->in1 = in1;
    obj->in2 = in2;
    obj->out = out;
}

void mod_classify_disconnect(mod_classify_obj *obj)
{

    obj->in1 = (msg_hops_obj *)NULL;
    obj->in2 = (msg_tracks_obj *)NULL;
    obj->out = (msg_categories_obj *)NULL;
}

void mod_classify_enable(mod_classify_obj *obj)
{

    obj->enabled = 1;
}

void mod_classify_disable(mod_classify_obj *obj)
{

    obj->enabled = 0;
}

mod_classify_cfg *mod_classify_cfg_construct(void)
{

    mod_classify_cfg *cfg;

    cfg = (mod_classify_cfg *)malloc(sizeof(mod_classify_cfg));

    return cfg;
}

void mod_classify_cfg_destroy(mod_classify_cfg *cfg)
{

    free((void *)cfg);
}

void mod_classify_cfg_printf(const mod_classify_cfg *cfg)
{
}
