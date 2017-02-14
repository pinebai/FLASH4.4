#if 0
    This file contains the indices of the various attributes used within the
    Pipeline unit.
#endif

#if 0
    Pipeline status

    PL_STATUS_SIZE : total number of pipeline status entries

    PL_STATUS_COMM : The strategy is to add +1 to the local processor slot, if
                     a send is being posted and to add -1 to the channel processor
                     slot, if a receive from this channel has completed. When
                     checking the global pipeline status, an all reduce with
                     summation operation on all processors is made. The result is
                     the number of active send/receive pairs on each processor.
                     If a zero vector is the result, we know that the overall
                     pipeline communication is done. However, there can still be
                     items sitting in the receive buffers and in the item buffer.

    PL_STATUS_RECV : For each processor, an entry of +1 means that at least one of
                     the receive buffers for the channels has still non-processed
                     items (items not stored away in the item buffer). An entry of
                     0 indicates thus that ALL receive buffers have been processed.

    PL_STATUS_ITEM : For each processor, an entry of +1 indicates that there are
                     still items left in the buffer, which have no been retrieved
                     by the external application. An entry of 0 means the item
                     buffer is empty, i.e. its item count is 0.

#endif

#define PL_STATUS_SIZE 3

#define PL_STATUS_COMM 1
#define PL_STATUS_RECV 2
#define PL_STATUS_ITEM 3

