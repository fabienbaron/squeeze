/*  Programme pour tester le generateur   RngStreams.c  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RngStream.h"

int main ()
{
   double sum = 0.0, sum3 = 0.0;
   long sumi = 0;
   int  i;
   RngStream g1, g2, g3;
   RngStream gar[4];
   unsigned long germe[6] = { 1, 1, 1, 1, 1, 1 };

   g1 = RngStream_CreateStream ("g1");
   g2 = RngStream_CreateStream ("g2");
   g3 = RngStream_CreateStream ("g3");

   printf ("Initial states of g1, g2, and g3:\n\n");
   RngStream_WriteState (g1);
   RngStream_WriteState (g2);
   RngStream_WriteState (g3);
   sum = RngStream_RandU01 (g2) + RngStream_RandU01 (g3);
   for (i = 0;  i < 12345; i++)
      RngStream_RandU01 (g2);

   RngStream_AdvanceState (g1, 5, 3);   
   printf ("State of g1 after advancing by 2^5 + 3 = 35 steps:\n");
   RngStream_WriteState (g1);
   printf ("RandU01 (g1) = %12.8f\n\n", RngStream_RandU01 (g1));

   RngStream_ResetStartStream (g1);
   for (i = 0;  i < 35; i++)    RngStream_AdvanceState (g1, 0,1);
   printf ( "State of g1 after reset and advancing 35 times by 1:\n");
   RngStream_WriteState (g1);
   printf ( "RandU01 (g1) = %12.8f\n\n", RngStream_RandU01 (g1));

   RngStream_ResetStartStream (g1);
   for (i = 0;  i < 35; i++)    sumi += RngStream_RandInt (g1, 1, 10);
   printf ( "State of g1 after reset and 35 calls to RandInt (1, 10):\n");
   RngStream_WriteState (g1);
   printf ( "   sum of 35 integers in [1, 10] = %ld\n\n", sumi);
   sum += sumi / 100.0;
   printf ( "RandU01 (g1) = %12.8f\n\n", RngStream_RandU01 (g1));

   sum3 = 0.0;
   RngStream_ResetStartStream (g1);
   RngStream_IncreasedPrecis (g1, 1);
   sumi = 0;
   for (i = 0;  i < 17; i++)     sumi += RngStream_RandInt (g1, 1, 10);
   printf ("State of g1 after reset, IncreasedPrecis (1) and 17 calls"
           " to RandInt (1, 10):\n");
   RngStream_WriteState (g1);
   RngStream_IncreasedPrecis (g1, 0);
   RngStream_RandInt (g1, 1, 10);
   printf ("State of g1 after IncreasedPrecis (0) and 1 call to RandInt\n");
   RngStream_WriteState (g1);
   sum3 = sumi / 10.0;

   RngStream_ResetStartStream (g1);
   RngStream_IncreasedPrecis (g1, 1);
   for (i = 0;  i < 17; i++)    sum3 += RngStream_RandU01 (g1);
   printf ("State of g1 after reset, IncreasedPrecis (1) and 17 calls"
           " to RandU01:\n");
   RngStream_WriteState (g1);
   RngStream_IncreasedPrecis (g1, 0);
   RngStream_RandU01 (g1);
   printf ("State of g1 after IncreasedPrecis (0) and 1 call to RandU01\n");
   RngStream_WriteState (g1);
   sum += sum3 / 10.0;

   sum3 = 0.0;
   printf ( "Sum of first 100 output values from stream g3:\n");
   for (i = 0;  i < 100;  i++) {
      sum3 += RngStream_RandU01 (g3);
   }
   printf ( "   sum = %12.6f\n\n", sum3);
   sum += sum3 / 10.0;

   printf ( "\nReset stream g3 to its initial seed.\n");
   RngStream_ResetStartSubstream (g3);
   printf ( "First 5 output values from stream g3:\n");
   for (i=1; i<=5; i++)
      printf ("%12.8f\n", RngStream_RandU01 (g3));
   sum += RngStream_RandU01 (g3);

   printf ("\nReset stream g3 to the next Substream, 4 times.\n");
   for (i=1; i<=4; i++)
      RngStream_ResetNextSubstream (g3);
   printf ("First 5 output values from stream g3, fourth Substream:\n");
   for (i=1; i<=5; i++)
      printf ("%12.8f\n", RngStream_RandU01 (g3));
   sum += RngStream_RandU01 (g3);

   printf ("\nReset stream g2 to the beginning of Substream.\n");
   RngStream_ResetStartSubstream (g2);
   printf (" Sum of 100000 values from stream g2 with double precision:   ");
   RngStream_IncreasedPrecis (g2, 1);
   sum3 = 0.0;
   for (i=1; i<=100000; i++)
      sum3 += RngStream_RandU01 (g2);
   printf ("%12.4f\n", sum3);
   sum += sum3 / 10000.0;
   RngStream_IncreasedPrecis (g2, 0);

   RngStream_SetAntithetic (g3, 1);
   printf (" Sum of 100000 antithetic output values from stream g3:   ");
   sum3 = 0.0;
   for (i=1; i<=100000; i++)
      sum3 += RngStream_RandU01 (g3);
   printf ("%12.4f\n", sum3);
   sum += sum3 / 10000.0;

   RngStream_DeleteStream (&g3);
   RngStream_DeleteStream (&g2);
   RngStream_DeleteStream (&g1);

   printf ("\nSetPackageSeed to seed = { 1, 1, 1, 1, 1, 1 }\n");
   RngStream_SetPackageSeed (germe);

   printf ("\nCreate an array of 4 named streams and "
           "write their full state\n");
   gar[0] = RngStream_CreateStream ("Poisson");
   gar[1] = RngStream_CreateStream ("Laplace");
   gar[2] = RngStream_CreateStream ("Galois");
   gar[3] = RngStream_CreateStream ("Cantor");

   for  (i = 0; i < 4; i++)
      RngStream_WriteStateFull (gar[i]);

   printf ("Jump stream Galois by 2^127 steps backward\n");
   RngStream_AdvanceState (gar[2], -127, 0);
   RngStream_WriteState (gar[2]);
   RngStream_ResetNextSubstream (gar[2]);

   for  (i = 0; i < 4; i++)
      sum += RngStream_RandU01 (gar[i]);

   printf ("--------------------------------------\n");
   printf ("Final Sum = %.6f", sum);

   if (fabs(sum-23.705324)<1e-6) {
     printf ("\t ... ok\n\n");
     exit (0);
   }
   else {
     printf ("\t ... failed\n\n");
     exit (1);
   }

}
