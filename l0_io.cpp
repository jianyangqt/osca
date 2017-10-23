//
//  l0_io.cpp
//  osc
//
//  Created by Futao Zhang on 31/03/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "l0_io.h"


void read_msglist(string msglistfile, vector<string> &msglist, string msg)
{
    // Read msglist file
    msglist.clear();
    string StrBuf;
    ifstream i_msglist(msglistfile.c_str());
    if(!i_msglist)
    {
        
        LOGPRINTF("Error: can not open the file %s to read.\n",msglistfile.c_str());
        TERMINATE();
    }
    LOGPRINTF( "Reading a list of %s from %s .\n", msg.c_str() ,msglistfile.c_str());
  
    while(i_msglist>>StrBuf){
        msglist.push_back(StrBuf);
        getline(i_msglist, StrBuf);
    }
    i_msglist.close();
}

void update_map_kp(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep)
{
    int i=0;
    map<string, int> id_map_buf(id_map);
    for(i=0; i<id_list.size(); i++) id_map_buf.erase(id_list[i]);
    map<string, int>::iterator iter;
    for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) id_map.erase(iter->first);
    
    keep.clear();
    for(iter=id_map.begin(); iter!=id_map.end(); iter++) keep.push_back(iter->second);
    stable_sort(keep.begin(), keep.end());
    /* 
     //the same as above
     keep.clear();
     map<string, int> id_map_buf;
     map<string, int>::iterator iter;
     for(int i=0;i<id_list.size();i++)
     {
     iter=id_map.find(id_list[i]);
     if(iter!=id_map.end()){
     keep.push_back(iter->second);
     id_map_buf.insert(pair<string,int>(id_list[i],iter->second));
     } else {
     cout<<id_list[i]<<endl;
     }
     }
     id_map=id_map_buf;
     stable_sort(keep.begin(), keep.end());
     printf("map: %ld, keep: %ld \n",id_map.size(),keep.size());
     */
}
void update_map_rm(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep)
{
    int i = 0;
    for (i = 0; i < id_list.size(); i++) id_map.erase(id_list[i]);
    
    keep.clear();
    map<string, int>::iterator iter;
    for (iter = id_map.begin(); iter != id_map.end(); iter++) keep.push_back(iter->second);
    stable_sort(keep.begin(), keep.end());
}
void read_indi_list(string indi_list_file, vector<string> &indi_list)
{
    ifstream i_indi_list(indi_list_file.c_str());
    if(!i_indi_list)
    {
        LOGPRINTF("ERROR: can not open the file %s to read.\n", indi_list_file.c_str());
        TERMINATE();
    }
    string str_buf, id_buf;
    indi_list.clear();
    while(i_indi_list){
        i_indi_list>>str_buf;
        if(i_indi_list.eof()) break;
        id_buf=str_buf+":";
        i_indi_list>>str_buf;
        id_buf+=str_buf;
        indi_list.push_back(id_buf);
        getline(i_indi_list, str_buf);
    }
    i_indi_list.close();
}
void write_msglist(char* outFileName, vector<string> &msg)
{
    FILE* efile=NULL;
    if(fopen_checked(&efile, outFileName,"w")) TERMINATE();
    
    for (int i = 0; i <msg.size(); i++) {
        string str=msg[i]+'\n';
        if(fputs_checked(str.c_str(),efile))
        {
            LOGPRINTF("ERROR: in writing file %s .\n", outFileName);
            TERMINATE();
        }
    }
    fclose(efile);
}
