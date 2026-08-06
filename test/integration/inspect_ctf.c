#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <samlib.h>

int main(int argc, char **argv)
{
    HeaderInfo header;
    ChannelInfo *channels;
    EpochInfo *epochs;
    SAM_MARKS *markers;
    unsigned char **bad;
    char dataset[1024];
    char marker_file[1200];
    int marker_count;
    int stim_count = 0;
    int missing_count = 0;
    int i;

    if (argc != 2) {
        fprintf(stderr, "usage: %s DATASET.ds\n", argv[0]);
        return EXIT_FAILURE;
    }
    if (strlen(argv[1]) >= sizeof(dataset)) {
        fprintf(stderr, "dataset path is too long\n");
        return EXIT_FAILURE;
    }
    strcpy(dataset, argv[1]);
    GetDsInfo(dataset, &header, &channels, &epochs, &bad, TRUE);

    snprintf(marker_file, sizeof(marker_file), "%s/MarkerFile.mrk", argv[1]);
    GetMarkers(marker_file, &markers, &marker_count);
    for (i = 0; i < marker_count; i++) {
        if (strcmp(markers[i].Name, "stim") == 0)
            stim_count++;
        else if (strcmp(markers[i].Name, "missingstim") == 0)
            missing_count++;
    }

    printf(
        "{\"epochs\":%d,\"primary_channels\":%d,\"samples\":%d,"
        "\"sample_rate\":%.9g,\"markers\":%d,\"stim\":%d,\"missingstim\":%d}\n",
        header.NumEpochs,
        header.NumPri,
        header.MaxSamples,
        header.SampleRate,
        marker_count,
        stim_count,
        missing_count
    );
    return EXIT_SUCCESS;
}
