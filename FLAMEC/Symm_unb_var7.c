/*

censing information see
http://www.cs.utexas.edu/users/flame/license.html 

Programmed by: Jay Shah, Matt Ebeweber
Email of author
*/

#include "FLAME.h"

int Symm_unb_var7( FLA_Obj A, FLA_Obj B, FLA_Obj C )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
          A20,  a21,     A22;

  FLA_Obj BT,              B0,
          BB,              b1t,
          B2;

  FLA_Obj CT,              C0,
          CB,              c1t,
          C2;

  FLA_Part_2x2( A,    &ATL, &ATR,
      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x1( B,    &BT, 
      &BB,            0, FLA_TOP );

  FLA_Part_2x1( C,    &CT, 
      &CB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
        /* ************* */   /* ************************** */
        &a10t, /**/ &alpha11, &a12t,
        ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
        1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
        /* ** */            /* *** */
        &b1t, 
        BB,                &B2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
        /* ** */            /* *** */
        &c1t, 
        CB,                &C2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    
    // C0 = C0 + a10*b1t;
	  FLA_Ger(FLA_ONE, a10t, b1t, C0); 

    // c1t = c1t + alpha11*b1t;
	  FLA_Axpy(alpha11, b1t, c1t); 

    // C2 = C2 + a21t*b1t;
    FLA_Ger(FLA_ONE, a21, b1t, C2);

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
        a10t, alpha11, /**/ a12t,
        /* ************** */  /* ************************ */
        &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
        FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
        b1t, 
        /* ** */           /* *** */
        &BB,                B2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
        c1t, 
        /* ** */           /* *** */
        &CB,                C2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}

