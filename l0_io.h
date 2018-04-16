//
//  l0_io.h
//  osc
//
//  Created by Futao Zhang on 31/03/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __osc__l0_io__
#define __osc__l0_io__

#include "l0_com.h"
#include "l0_mem.h"

extern void read_msglist(string msglistfile, vector<string> &msglist, string msg);
extern void update_map_kp(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep);
extern void update_map_rm(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep);
extern void read_indi_list(string indi_list_file, vector<string> &indi_list);
extern void write_msglist(char* outFileName, vector<string> &msg);
extern uint32_t readuint32(FILE *f);
uint64_t readuint64(FILE *f);
float readfloat(FILE *f);
#endif /* defined(__osc__l0_io__) */
