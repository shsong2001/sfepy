#include "terms.h"
#include "geommech.h"

#undef __FUNC__
#define __FUNC__ "dw_surface_v_dot_n_s"
int32 dw_surface_v_dot_n_s(FMField *out,
                           FMField *coef, FMField *val_qp,
                           Mapping *rsg,
                           Mapping *csg,
                           int32 isDiff)
{
  int32 ii, dim, nQP, nEPR, nEPC, ret = RET_OK;
  FMField *aux1 = 0, *aux2 = 0;

  nQP = rsg->normal->nLev;
  dim = rsg->normal->nRow;
  nEPR = rsg->bf->nCol;
  nEPC = csg->bf->nCol;

  fmf_createAlloc(&aux1, 1, nQP, dim * nEPR, 1);
  if (isDiff) {
    fmf_createAlloc(&aux2, 1, nQP, dim * nEPR, nEPC);
  } else {
    fmf_createAlloc(&aux2, 1, nQP, dim * nEPR, 1);
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell(out, ii);
    FMF_SetCellX1(coef, ii);
    FMF_SetCell(rsg->det, ii);
    FMF_SetCellX1(rsg->bf, ii);

    if (isDiff) {
      FMF_SetCellX1(csg->bf, ii);
      bf_actt(aux1, rsg->bf, csg->normal);
      fmf_mulAB_nn(aux2, aux1, csg->bf);
      fmf_mul(aux2, coef->val);
      fmf_sumLevelsMulF(out, aux2, rsg->det->val);
    } else {
      FMF_SetCell(val_qp, ii);
      bf_actt(aux1, rsg->bf, csg->normal);
      fmf_mulAB_nn(aux2, aux1, val_qp);
      fmf_mul(aux2, coef->val);
      fmf_sumLevelsMulF(out, aux2, rsg->det->val);
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy(&aux1);
  fmf_freeDestroy(&aux2);

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_surface_s_v_dot_n"
int32 dw_surface_s_v_dot_n(FMField *out,
                           FMField *coef, FMField *val_qp,
                           Mapping *rsg,
                           Mapping *csg,
                           int32 isDiff)
{
  int32 ii, dim, nQP, nEPR, nEPC, ret = RET_OK;
  FMField *aux1 = 0, *aux2 = 0;

  nQP = rsg->det->nLev;
  dim = csg->normal->nRow;
  nEPR = rsg->bf->nCol;
  nEPC = csg->bf->nCol;

  if (isDiff) {
    fmf_createAlloc(&aux2, 1, nQP, nEPR, dim * nEPC);
    fmf_createAlloc(&aux1, 1, nQP, dim * nEPC, 1);
  } else {
    fmf_createAlloc(&aux2, 1, nQP, nEPR, 1);
    fmf_createAlloc(&aux1, 1, nQP, 1, 1);
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell(out, ii);
    FMF_SetCellX1(coef, ii);
    FMF_SetCell(rsg->det, ii);
    FMF_SetCellX1(rsg->bf, ii);

    if (isDiff) {
      FMF_SetCellX1(csg->bf, ii);
      bf_actt(aux1, csg->bf, csg->normal);
      fmf_mulATBT_nn(aux2, rsg->bf, aux1);
      fmf_mul(aux2, coef->val);
      fmf_sumLevelsMulF(out, aux2, rsg->det->val);
    } else {
      FMF_SetCell(val_qp, ii);
      fmf_mulATB_nn(aux1, csg->normal, val_qp);
      fmf_mulATB_nn(aux2, rsg->bf, aux1);
      fmf_mul(aux2, coef->val);
      fmf_sumLevelsMulF(out, aux2, rsg->det->val);
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy(&aux1);
  fmf_freeDestroy(&aux2);

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_volume_dot_vector"
/*!
  @par Revision history:
  - 21.11.2006, c
*/
int32 dw_volume_dot_vector( FMField *out, FMField *coef, FMField *val_qp,
                            Mapping *rvg, Mapping *cvg,
                            int32 isDiff )
{
  printf("========================================\n");
  printf("C JV: terms_dot.c : dw_volume_dot_vector\n");
  int32 ii, dim, nc, nQP, nEPR, nEPC, ret = RET_OK;
  FMField *ftfu = 0, *ftf1 = 0, *ftf = 0, *cf = 0, *cfu = 0;

  nQP = rvg->nQP;
  dim = rvg->dim;
  nEPR = rvg->bf->nCol;
  nEPC = cvg->bf->nCol;
  nc = coef->nCol;
  /*printf("coef->(nCell, nLev, nRow, nCol), (nAlloc, cellSize)\n");
  printf("(%d, %d, %d, %d), (%d, %d)\n", coef->nCell, coef->nLev, coef->nRow,
    		 coef->nCol, coef->nAlloc, coef->cellSize);
  printf("rvg->bf->(nCell, nLev, nRow, nCol), (nAlloc, cellSize)\n");
  printf("(%d, %d, %d, %d), (%d, %d)\n", rvg->bf->nCell, rvg->bf->nLev,
		  rvg->bf->nRow, rvg->bf->nCol, rvg->bf->nAlloc, rvg->bf->cellSize);
  printf("cvg->bf->(nCell, nLev, nRow, nCol), (nAlloc, cellSize)\n");
  printf("(%d, %d, %d, %d), (%d, %d)\n", cvg->bf->nCell, cvg->bf->nLev,
		  cvg->bf->nRow, cvg->bf->nCol, cvg->bf->nAlloc, cvg->bf->cellSize);
  printf("out->(nCell, nLev, nRow, nCol), (nAlloc, cellSize)\n");
  printf("(%d, %d, %d, %d), (%d, %d)\n", out->nCell, out->nLev, out->nRow,
    		 out->nCol, out->nAlloc, out->cellSize);*/

  int32 nRow, nCol, switcher, Col_switcher, Row_switcher;
  Row_switcher = 0;
  Col_switcher = 0;
  if (rvg->bf->nRow > 1) {
	  nRow = nEPR;
	  Row_switcher = 1;
  } else{
	  nRow = nEPR * dim;
  }
  if (cvg->bf->nRow > 1) {
  	  nCol = nEPC;
  	  Col_switcher = 1;
    } else{
  	  nCol = nEPC * dim;
    }
  switcher = Row_switcher + Col_switcher;

  /*printf("nRow=%d; nCol=%d; switcher=%d; \n", nRow, nCol, switcher);
  printf("isDiff=%d; nc=%d; \n", isDiff, nc);*/

  if (isDiff) {
    /*fmf_createAlloc( &ftf, 1, nQP, nEPR * dim, nEPC * dim );*/
    fmf_createAlloc( &ftf, 1, nQP, nRow, nCol );
    if (nc == 1) {
      fmf_createAlloc( &ftf1, 1, nQP, nEPR, nEPC );
    } else {
      /*fmf_createAlloc( &cf, 1, nQP, dim, dim * nEPC );*/
      fmf_createAlloc( &cf, 1, nQP, dim, nRow );
    }
  } else {
    /*fmf_createAlloc( &ftfu, 1, nQP, dim * nEPR, 1 );*/
    fmf_createAlloc( &ftfu, 1, nQP, nRow, 1 );
    if (nc > 1) {
      fmf_createAlloc( &cfu, 1, nQP, dim, 1 );
    }
  }

  int32 iqp, ic, ir, kk, jj, nRC, nCd;
  float64 mat, u, v;
  nRC = nCol*nRow;
  nCd = nCol*dim;
  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCellX1( coef, ii );
    FMF_SetCell( rvg->det, ii );
    FMF_SetCellX1( rvg->bf, ii );

    if (isDiff) {
      FMF_SetCellX1( cvg->bf, ii );

      if (nc == 1) {
    	if (switcher == 0){
        fmf_mulATB_nn( ftf1, rvg->bf, cvg->bf );
        bf_buildFTF( ftf, ftf1 );
        fmf_mul( ftf, coef->val );
    	}else if (Col_switcher == 0){
    	  /*printf("term 03 \n");
    	  printf("rvg->bf:\n");
    	  fmf_print( rvg->bf, stdout, 0 );
    	  printf("nQP=%d, nCol=%d, nRow=%d, nEPR, nRC \n",
    			  nQP, nCol, nRow, nEPR, nRC);*/
    	  fmf_fillC( ftf, 0.0 );
    	  for (iqp=0; iqp<nQP; iqp++){
    		mat = coef->val[iqp];
    		  for (ic=0; ic<nCol; ic++){
  			    for (ir=0; ir<nRow; ir++){
    		      u = rvg->bf->val[iqp*nRC + ic*nRow + ir];
    		      ftf->val[iqp*nRC + ir*nCol + ic] += mat*u;
    		      /*printf("indices: iqp=%d, ic=%d, ir=%d, iu=%d, iftf=%d \n",
    		    		  iqp, ic, ir, iqp*nCd + ic*nCol + ir,
    		    		  iqp*nRC + ir*nCol + ic);
    		      printf("mat=%f, u=%f, res=%f \n", mat, u, mat*u);
    		      printf("ftf:\n");
    		      fmf_print( ftf, stdout, 0 );*/
    		  }
    		}
    	  }
    	}
    	else if (Row_switcher == 0){
    	  /*printf("term 30 \n");
    	  printf("cvg->bf:\n");
    	  fmf_print( cvg->bf, stdout, 0 );
    	  printf("nQP=%d, nCol=%d, nRow=%d, nEPR, nRC \n",
    	      	  nQP, nCol, nRow, nEPR, nRC);*/
      	  fmf_fillC( ftf, 0.0 );
      	  for (iqp=0; iqp<nQP; iqp++){
      		mat = coef->val[iqp];
      		for (ic=0; ic<nCol; ic++){
    		  for (ir=0; ir<nRow; ir++){
      		    u = cvg->bf->val[iqp*nRC + ir*nCol + ic];
      		    ftf->val[iqp*nRC + ir*nCol + ic] += mat*u;
      		    /*printf("indices: iqp=%d, ic=%d, ir=%d, iu=%d, iftf=%d \n",
      		    		  iqp, ic, ir, iqp*nCd + ic*nCol + ir,
      		    		  iqp*nRC + ir*nCol + ic);
      		    printf("mat=%f, u=%f, res=%f \n", mat, u, mat*u);
      		    printf("ftf:\n");
      		    fmf_print( ftf, stdout, 0 );*/
      		  }
      		}
      	  }
    	}
      } else {
    	if (switcher == 0){
        bf_ract( cf, cvg->bf, coef );
        bf_actt( ftf, rvg->bf, cf );
    	}else if (switcher == 2){
        fmf_fillC( ftf, 0.0 );
        for (iqp=0; iqp<nQP; iqp++){
          for (ir=0; ir<nRow; ir++){
        	for (ic=0; ic<nCol; ic++){
        	  for (kk=0; kk<dim; kk++){
        		for (jj=0; jj<dim; jj++){
        		  mat = coef->val[iqp*dim*dim + dim*kk + jj];
        		  u = rvg->bf->val[iqp*nCd + kk*nCol + ir];
        		  v = cvg->bf->val[iqp*nCd + jj*nCol + ic];
        		  ftf->val[iqp*nRC + ir*nCol + ic] += mat*u*v;
        		}
        	  }
        	}
          }
        }
    	}
      }
      fmf_sumLevelsMulF( out, ftf, rvg->det->val );
    } else {
      FMF_SetCell( val_qp, ii );

      if (nc == 1) {
        bf_actt( ftfu, rvg->bf, val_qp );
        fmf_mul( ftfu, coef->val );
      } else {
        fmf_mulAB_nn( cfu, coef, val_qp );
        bf_actt( ftfu, rvg->bf, cfu );
      }

      fmf_sumLevelsMulF( out, ftfu, rvg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftf );
    if (nc == 1) {
      fmf_freeDestroy( &ftf1 );
    } else {
      fmf_freeDestroy( &cf );
    }
  } else {
    fmf_freeDestroy( &ftfu );
    if (nc > 1) {
      fmf_freeDestroy( &cfu );
    }
  }
  /*printf("========================================\n");*/
  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_volume_dot_scalar"
/*!
  @par Revision history:
  - 01.02.2008, c
*/
int32 dw_volume_dot_scalar( FMField *out, FMField *coef, FMField *val_qp,
                            Mapping *rvg, Mapping *cvg,
                            int32 isDiff )
{
  int32 ii, nQP, nEPR, nEPC, ret = RET_OK;
  FMField *ftfp = 0, *ftf = 0, *cftf = 0;

  nQP = rvg->nQP;
  nEPR = rvg->bf->nCol;
  nEPC = cvg->bf->nCol;

  if (isDiff) {
    fmf_createAlloc( &ftf, 1, nQP, nEPR, nEPC );
    fmf_createAlloc( &cftf, 1, nQP, nEPR, nEPC );
  } else {
    fmf_createAlloc( &ftfp, 1, nQP, nEPR, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( rvg->det, ii );
    FMF_SetCellX1( coef, ii );
    FMF_SetCellX1( rvg->bf, ii );

    if (isDiff) {
      FMF_SetCellX1( cvg->bf, ii );

      fmf_mulATB_nn( ftf, rvg->bf, cvg->bf );
      fmf_mulAF( cftf, ftf, coef->val );
      fmf_sumLevelsMulF( out, cftf, rvg->det->val );
    } else {
      FMF_SetCell( val_qp, ii );

      bf_actt( ftfp, rvg->bf, val_qp );
      fmf_mul( ftfp, coef->val );
      fmf_sumLevelsMulF( out, ftfp, rvg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftf );
    fmf_freeDestroy( &cftf );
  } else {
    fmf_freeDestroy( &ftfp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_surface_dot_scalar"
/*!
  @par Revision history:
  - 09.03.2009, c
*/
#undef __FUNC__
#define __FUNC__ "dw_v_dot_grad_s_vw"
int32 dw_v_dot_grad_s_vw( FMField *out, FMField *coef, FMField *grad,
                          Mapping *vvg, Mapping *svg,
                          int32 isDiff )
{
  int32 ii, nc, nEPV, nEPS, dim, nQP, ret = RET_OK;
  FMField *ftg = 0, *cg = 0;

  nQP = vvg->bfGM->nLev;
  dim = vvg->bfGM->nRow;
  nEPS = svg->bfGM->nCol;
  nEPV = vvg->bf->nCol;
  nc = coef->nCol;

  if (isDiff == 1) {
    fmf_createAlloc( &ftg, 1, nQP, dim * nEPV, nEPS );
    if (nc > 1) {
      fmf_createAlloc( &cg, 1, nQP, dim, nEPS );
    }
  } else {
    fmf_createAlloc( &ftg, 1, nQP, dim * nEPV, 1 );
    if (nc > 1) {
      fmf_createAlloc( &cg, 1, nQP, dim, 1 );
    }
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCellX1( coef, ii );
    FMF_SetCell( vvg->det, ii );
    FMF_SetCellX1( vvg->bf, ii );

    if (isDiff == 1) {
      FMF_SetCell( svg->bfGM, ii );

      if (nc == 1) {
        bf_actt( ftg, vvg->bf, svg->bfGM );
        fmf_mul( ftg, coef->val );
      } else {
        // Phi^T C Gc
        fmf_mulAB_nn( cg, coef, svg->bfGM );
        bf_actt( ftg, vvg->bf, cg );
      }
    } else {
      FMF_SetCell( grad, ii );

      if (nc == 1) {
        bf_actt_c1( ftg, vvg->bf, grad );
        fmf_mul( ftg, coef->val );
      } else {
        // Phi^T C Gc s
        fmf_mulAB_nn( cg, coef, grad );
        bf_actt( ftg, vvg->bf, cg );
      }
    }
    fmf_sumLevelsMulF( out, ftg, vvg->det->val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &ftg );
  fmf_freeDestroy( &cg );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_v_dot_grad_s_sw"
int32 dw_v_dot_grad_s_sw( FMField *out, FMField *coef, FMField *val_qp,
                          Mapping *vvg, Mapping *svg,
                          int32 isDiff )
{
  int32 ii, nc, nEPV, nEPS, dim, nQP, ret = RET_OK;
  FMField *gtf = 0, *ctf = 0, *ct = 0;

  nQP = vvg->bfGM->nLev;
  dim = vvg->bfGM->nRow;
  nEPS = svg->bfGM->nCol;
  nEPV = vvg->bf->nCol;
  nc = coef->nCol;

  if (isDiff == 1) {
    fmf_createAlloc( &gtf, 1, nQP, nEPS, dim * nEPV );
    if (nc > 1) {
      fmf_createAlloc( &ctf, 1, nQP, dim, dim * nEPV );
      fmf_createAlloc( &ct, 1, nQP, dim, dim );
    } else {
      // Gc^T.
      fmf_createAlloc( &ctf, 1, nQP, nEPS, dim );
    }
  } else {
    fmf_createAlloc( &gtf, 1, nQP, nEPS, 1 );
    if (nc > 1) {
      fmf_createAlloc( &ctf, 1, nQP, dim, 1 );
    }
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCellX1( coef, ii );
    FMF_SetCell( svg->bfGM, ii );
    FMF_SetCell( vvg->det, ii );

    if (isDiff == 1) {
      FMF_SetCellX1( vvg->bf, ii );

      if (nc == 1) {
        // Transpose Gc.
        fmf_mulATC( ctf, svg->bfGM, 1.0 );
        // Gc^T Phi.
        bf_ract( gtf, vvg->bf, ctf );
        fmf_mul( gtf, coef->val );
      } else {
        // Transpose C.
        fmf_mulATC( ct, coef, 1.0 );
        // Gc^T C^T Phi.
        bf_ract( ctf, vvg->bf, ct );
        fmf_mulATB_nn( gtf, svg->bfGM, ctf );
      }
    } else {
      FMF_SetCell( val_qp, ii );

      if (nc == 1) {
        fmf_mulATB_nn( gtf, svg->bfGM, val_qp );
        fmf_mul( gtf, coef->val );
      } else {
        // Gc^T C^T Phi v.
        fmf_mulATB_nn( ctf, coef, val_qp );
        fmf_mulATB_nn( gtf, svg->bfGM, ctf );
      }
    }
    fmf_sumLevelsMulF( out, gtf, vvg->det->val );
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &gtf );
  fmf_freeDestroy( &ctf );
  fmf_freeDestroy( &ct );

  return( ret );
}
