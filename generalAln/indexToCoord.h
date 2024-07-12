/*########################################################
# Name indexToCoord
# Use:
#  o Holds functions to convert alignment matrix index's
#    to coordinates
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'   o header:
'     - Included libraries
'   o fun01: indexToQry_water
'     - Gets the query coordinates of the query sequence
'       in an matrix.
'   o fun02: indexToRef_water
'     - Gets the coordinates of the reference sequence in
'       in an matrix.
'   o fun03: indexToCoord_water
'     - Gets the coordinates of the reference and query
'       sequence in an matrix.
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards
\-------------------------------------------------------*/

#ifndef INDEX_TO_COORDINATES_H
#define INDEX_TO_COORDINATES_H

/*-------------------------------------------------------\
| Fun01: indexToQry_indexToCoord
|   - Gets the query coordinates of the query sequence in
|     an matrix.
| Input
|   - refLen:
|     o Length of the reference sequence
|   - index:
|     o Index to convert to cooridnates
|   - qryCoord:
|     o Will hold the coordinate of the query sequence
| Output:
|  - Returns:
|    o The query coordinate that was in the index
\-------------------------------------------------------*/
#define \
indexToQry_indexToCoord( \
   refLen, \
   index \
)({ \
   ulong resMacUL = 0; \
   resMacUL = ( (index) / ((refLen) + 1) ); \
   resMacUL += (resMacUL == 0); \
   resMacUL -= ( ( (index) % ((refLen) + 1) ) > 0 ); \
   /* Find query coordinates from an index:
   `  - refLen + 1:
   `    o gives the length of each row in the matrix
   `    o The + 1 accounts for the gap column
   `  - index / (refLen + 1)
   `    o The number of rows down, which is the number of
   `      query bases + the gap row
   `  - index >= refLen
   `    o Is 1 if the index is divisible by the reference
   `      This means I have an index 1 value
   `    o else index is in the gap row (index 0)
   `  - position - (index 1 or 0)
   `    o Removes the gap row from the query if it is
   `      present. This makes sure I get an index 0 result
   */\
}) /*indexToQry_indexToCoord*/

/*-------------------------------------------------------\
| Fun02: indexToRef_indexToCoord
|   - Gets the coordinates of the reference sequence in
|     an matrix.
| Input
|   - refLen:
|     o Length of the reference sequence
|   - index:
|     o Index to convert to cooridnates
|   - refCoord:
|     o Will hold the coordinate of the reference sequence
| Output:
|  - Returns
|    o The reference coordinate that was in the index
\-------------------------------------------------------*/
#define \
indexToRef_indexToCoord( \
   refLen, \
   index \
)({ \
   ulong resMacUL = 0; \
   resMacUL = ( (index) % ((refLen) + 1) ); \
   resMacUL += (resMacUL == 0); \
   resMacUL -= ( ( (index) / ((refLen) + 1) ) > 0 ); \
   resMacUL; \
   /* Find reference coordinates from an index:
   `  - refLen + 1:
   `    o gives the length of each row in the matrix
   `    o The + 1 accounts for the gap column
   `  - index % (refLen + 1)
   `    o The number of columns, which is the number of
   `      reference bases + the gap row
   `  - index / (refLen + 1)
   `    o Gets the query position
   `  - (index / (refLen + 1)) > 0
   `    o Is 0 if I am in the gap column, 1 if on a real
   `      reference base
   `  - position - (index 1 or 0)
   `    o Removes the gap column from the reference if it
   `      is present. This makes sure I get an index 0
   `      result.
   */\
}) /*indexToRef_indexToCoord*/

/*-------------------------------------------------------\
| Fun03: indexToCoord_indexToCoord
|   - Gets the coordinates of the reference and query
|     sequence in an matrix.
| Input
|   - refLen:
|     o Length of the reference sequence
|   - index:
|     o Index to convert to cooridnates
|   - refCoord:
|     o Will hold the coordinate of the reference sequence
|   - qryCoord:
|     o Will hold the coordinate of the query sequence
| Output:
|  - Sets
|    o refCoord to the reference coordinate in index
|    o qryCoord to the query coordinate in index
\-------------------------------------------------------*/
#define \
indexToCoord( \
   refLen, \
   index, \
   refCoord, \
   qryCoord \
){ \
   ulong resMacUL = 0; \
   \
   resMacUL = ( (index) % ((refLen) + 1) ); \
   (qryCoord) = ( (index) / ((refLen) + 1) ); \
   \
   (refCoord) = resMacUL; \
   (refCoord) += ((refCoord) == 0); \
   (refCoord) -= ((qryCoord) > 0); \
   \
   (qryCoord) += ((qryCoord) == 0); \
   (qryCoord) -= (resMacUL > 0); \
} /*indexToCoord_indexToCoord*/

#endif
