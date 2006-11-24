/******************************************************************/

#ifndef PA_MACROS_H
#define PA_MACROS_H

/******************************************************************/
// macro's for function derivs - variables
#define S y[0]
#define I y[1]
#define R y[2]
#define s y[3]
#define i y[4]
#define r y[5]
//---------------- Virus net
#define SS_d y[6]
#define SI_d y[7]
#define IS_d y[7]
#define SR_d y[8]
#define II_d y[9]
#define IR_d y[10]
#define RI_d y[10]
#define RR_d y[11]
//----------------
#define ss_d y[12]
#define si_d y[13]
#define is_d y[13]
#define sr_d y[14]
#define rs_d y[14]
#define ii_d y[15]
#define ir_d y[16]
#define ri_d y[16]
#define rr_d y[17]
//----------------
#define Ss_d y[18]
#define sS_d y[18]
#define Si_d y[19]
#define iS_d y[19]
#define Sr_d y[20]
#define rS_d y[20]
#define Ii_d y[21]
#define iI_d y[21]
#define Ir_d y[22]
#define rI_d y[22]
#define Rr_d y[23]
#define rR_d y[23]
//----------------
#define sI_d y[24]
#define Is_d y[24]
#define sR_d y[25]
#define Rs_d y[25]
#define iR_d y[26]
#define Ri_d y[26]

//---------------- Info net
#define SS_i y[27]
#define SI_i y[28]
#define IS_i y[28]
#define SR_i y[29]
#define II_i y[30]
#define IR_i y[31]
#define RI_i y[31]
#define RR_i y[32]
//----------------
#define ss_i y[33]
#define si_i y[34]
#define is_i y[34]
#define sr_i y[35]
#define rs_i y[35]
#define ii_i y[36]
#define ir_i y[37]
#define ri_i y[37]
#define rr_i y[38]
//----------------
#define Ss_i y[39]
#define sS_i y[39]
#define Si_i y[40]
#define iS_i y[40]
#define Sr_i y[41]
#define rS_i y[41]
#define Ii_i y[42]
#define iI_i y[42]
#define Ir_i y[43]
#define rI_i y[43]
#define Rr_i y[44]
#define rR_i y[44]
//----------------
#define sI_i y[45]
#define Is_i y[45]
#define sR_i y[46]
#define Rs_i y[46]
#define iR_i y[47]
#define Ri_i y[47]

// macro's for function derivs - d/dt
#define S_t rhs[0]
#define I_t rhs[1]
#define R_t rhs[2]
#define s_t rhs[3]
#define i_t rhs[4]
#define r_t rhs[5]
//---------------- Virus net
#define SS_d_t rhs[6]
#define SI_d_t rhs[7]
#define SR_d_t rhs[8]
#define II_d_t rhs[9]
#define IR_d_t rhs[10]
#define RR_d_t rhs[11]
//----------------
#define ss_d_t rhs[12]
#define si_d_t rhs[13]
#define sr_d_t rhs[14]
#define ii_d_t rhs[15]
#define ir_d_t rhs[16]
#define rr_d_t rhs[17]
//----------------
#define Ss_d_t rhs[18]
#define Si_d_t rhs[19]
#define Sr_d_t rhs[20]
#define Ii_d_t rhs[21]
#define Ir_d_t rhs[22]
#define Rr_d_t rhs[23]
#define rR_d_t rhs[23]
//----------------
#define sI_d_t rhs[24]
#define sR_d_t rhs[25]
#define iR_d_t rhs[26]

//---------------- Info net
#define SS_i_t rhs[27]
#define SI_i_t rhs[28]
#define SR_i_t rhs[29]
#define II_i_t rhs[30]
#define IR_i_t rhs[31]
#define RR_i_t rhs[32]
//----------------
#define ss_i_t rhs[33]
#define si_i_t rhs[34]
#define sr_i_t rhs[35]
#define ii_i_t rhs[36]
#define ir_i_t rhs[37]
#define rr_i_t rhs[38]
//----------------
#define Ss_i_t rhs[39]
#define Si_i_t rhs[40]
#define Sr_i_t rhs[41]
#define Ii_i_t rhs[42]
#define Ir_i_t rhs[43]
#define Rr_i_t rhs[44]
#define rR_i_t rhs[44]
//----------------
#define sI_i_t rhs[45]
#define sR_i_t rhs[46]
#define iR_i_t rhs[47]

/******************************************************************/

#endif // PA_MACROS_H
