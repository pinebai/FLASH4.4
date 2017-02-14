#include "io_ncmpi_name_to_id.h"

void FTOC(io_ncmpi_name_to_id)(const int * const pFileID,
			       const char name[],
			       const int * const pLen,
			       int *pVarID)
{
  const int fileID = *pFileID;
  const int len = *pLen;
  char name_c[MAX_STRING_LENGTH+1];

  assert(len > 0 && len <= MAX_STRING_LENGTH);
  strncpy(name_c, name, len);
  name_c[len] = '\0';

  io_ncmpi_name_to_id_c(fileID, name_c, pVarID);
}


void io_ncmpi_name_to_id_c(const int fileID,
			   const char name[],
			   int *pVarID)
{
  size_t nameLen;
  int err;
  char legacyName[LEGACY_LABEL_LEN+1];
  const char *underscorePad = "____"; /* 4 characters */

  err = ncmpi_inq_varid(fileID, name, pVarID);
  /* If the inquiry fails then we make another inquiry with a
     backwards-compatible underscore padded name.  We do this because
     the original I/O in FLASH creates variables which have names of 4
     characters in length, even if the underscore trimmed name is less
     than 4 characters, e.g. 'al' is 'al__'  */
  if (err == NC_ENOTVAR) {
    nameLen = strlen(name);
    if (nameLen < LEGACY_LABEL_LEN) {
      /* The if statement above ensures that the strncpy call below
	 will copy a null character into legacyName */
      strncpy(legacyName, name, LEGACY_LABEL_LEN);
      strncat(legacyName, underscorePad, LEGACY_LABEL_LEN - nameLen);
      err = ncmpi_inq_varid(fileID, legacyName, pVarID);
    }
  }
  assert(err == NC_NOERR);
}
