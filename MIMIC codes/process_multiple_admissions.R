rm(list = ls()) 

load(file="multiple_admissions.RData")
CHARTEVENTS <- read.csv("1000subjects/CHARTEVENTS10000.csv",header = T)
D_ITEMS <- read.csv("/home/ICU_data/mimic-iii-clinical-database-1.4/D_ITEMS.csv",header = T)
ICUSTAYS <- read.csv("1000subjects/ICUSTAYS10000.csv",header = T)
# CHARTEVENTS2 <- CHARTEVENTS[is.na(CHARTEVENTS$VALUE)==FALSE&is.na(CHARTEVENTS$VALUENUM)==FALSE,]
# heart rate: 211; Arterial Blood Pressure systolic,220050; Arterial Blood Pressure diastolic, 220051;
# respiratory rate, 220210; Arterial O2 Saturation, 220227; 

subject_ID<-sort(unique(CHARTEVENTS$SUBJECT_ID))
all_items<-unique(CHARTEVENTS$ITEMID)
filtered_subject_ID <- subject_ID
filtered_items<-all_items
# 
# i<-0
# for (item in all_items)   {
#   tmp_s<-unique(CHARTEVENTS2[CHARTEVENTS2$itemid==item,]$subject_id)
#   if (length(tmp_s)>=30)  {
#     i<-i+1
#     filtered_items[i]<-item
#   }
#   # print(c('subject:',as.character(i)))
#   # print(length(tmp_s))
# }
# filtered_items<-filtered_items[1:i]
# 
# 
# current_items<-all_items
# i<-0
# for (s in subject_ID)   {
#   tmp_item<-CHARTEVENTS2[CHARTEVENTS2$subject_id==s,]$itemid
#   if (length(unique(tmp_item))>=100)  {
#     i<-i+1
#     filtered_subject_ID[i]<-s
#   }
#   # current_items<-intersect(current_items,tmp_item)
#   # # print(i)
#   # print(length(unique(tmp_item)))
# }
# filtered_subject_ID<-filtered_subject_ID[1:i]
# 
# # partial subject records, items num >=100, subjects num>=30
# CHARTEVENTS3<-CHARTEVENTS2[(CHARTEVENTS2$subject_id%in%filtered_subject_ID)&(CHARTEVENTS2$itemid%in%filtered_items),]
# CHARTEVENTS3<-CHARTEVENTS3[is.na(CHARTEVENTS3$icustay_id)==FALSE,]
# CHARTEVENTS3<-cbind(CHARTEVENTS3,CHARTEVENTS3$stopped)
# colnames(CHARTEVENTS3)[16]<-'chartHour'
# for (i in 1:nrow(CHARTEVENTS3))   {
#   CHARTEVENTS3[i,]$chartHour<-difftime(CHARTEVENTS3[i,]$charttime,
#                                        ICUSTAYS[ICUSTAYS$icustay_id==CHARTEVENTS3[i,]$icustay_id,]$intime,
#                                        units="hours")
#   
# }
# 


CHARTEVENTS4<-CHARTEVENTS[is.na(CHARTEVENTS$ICUSTAY_ID)==FALSE,]    # all subjects
CHARTEVENTS4<-cbind(CHARTEVENTS4,CHARTEVENTS4$STOPPED)
colnames(CHARTEVENTS4)[16]<-'chartHour'
t1<-Sys.time()
for (i in 1:nrow(CHARTEVENTS4))   {
  CHARTEVENTS4[i,]$chartHour<-difftime(CHARTEVENTS4[i,]$CHARTTIME,
                                       ICUSTAYS[ICUSTAYS$ICUSTAY_ID==CHARTEVENTS4[i,]$ICUSTAY_ID,]$INTIME,
                                       units="hours")
  if(i%%50==0) {
    t2<-Sys.time()
    deltat<-difftime(t2,t1,units="secs")
    print(paste(as.character(i),', ',as.character(deltat)))
  } 
  
}


CHARTEVENTS4 <- read.csv("1000subjects/CHARTEVENTS10000_4.csv",header = T)
CHARTEVENTS4[CHARTEVENTS4$chartHour<0,]$chartHour<-0
subject_IDs<- unique(CHARTEVENTS4$SUBJECT_ID)


sj<-0
sw<-0
for (s in subject_IDs)   {
  sj<-sj+1
  print(sj)
  tmp_CHARTEVENTS<-CHARTEVENTS4[CHARTEVENTS4$SUBJECT_ID==s,]
  max_hour<-floor(as.numeric(max(tmp_CHARTEVENTS$chartHour)))
  min_hour<-abs(floor(as.numeric(min(tmp_CHARTEVENTS$chartHour))))
  #print(max_hour-min_hour)
  sub_item<-unique(tmp_CHARTEVENTS$ITEMID)
  s_matrix<-array(data=NA,dim=c(length(sub_item),max_hour-min_hour+1))
  if (nrow((s_matrix))<10 |ncol((s_matrix))<10)  next
  rownames(s_matrix)<-sub_item
  colnames(s_matrix)<-c(min_hour:max_hour)
  s_filtered_matrix <- s_matrix
  j <- 0 
  for (i in 1:nrow(s_matrix))  {
    for (t in min_hour:max_hour)  {
      tmp_record <- tmp_CHARTEVENTS[tmp_CHARTEVENTS$ITEMID==(rownames(s_matrix)[i])&
                                      abs(floor(as.numeric(tmp_CHARTEVENTS$chartHour)))==t,]
      if (nrow(tmp_record)>0)   {
        s_matrix[i,t-min_hour+1] <- mean(as.numeric(tmp_record$VALUE))
        #print(mean(as.numeric(tmp_record$value)))
      }
      #s_matrix[rownames(s_matrix)==tmp_CHARTEVENTS[i,]$itemid,floor(as.numeric(tmp_CHARTEVENTS[i,]$chartHour))+1]<-tmp_CHARTEVENTS[i,]$valuenum
    }
    if(length(s_matrix[i,is.na(s_matrix[i,])==FALSE]) >= ncol(s_matrix)*0.5 )   {
      j<- j+1
      s_filtered_matrix[j,]<-s_matrix[i,]
    }
  }
  if(j<=5) next
  s_filtered_matrix <- s_filtered_matrix[1:j,]
  if (nrow((s_filtered_matrix))<10 |ncol((s_filtered_matrix))<10)  next
  
  s_filtered_matrix2 <- s_filtered_matrix
  for (i in 1:nrow(s_filtered_matrix2))   {
    for (j in 1:3)    {
      nonNArecords <- s_filtered_matrix2[i,is.na(s_filtered_matrix2[i,])==FALSE]
      if (is.na(s_filtered_matrix2[i,j])) s_filtered_matrix2[i,j] <- mean(nonNArecords[1:3])
    }
    for (j in 4:ncol(s_filtered_matrix2))   {
      if (is.na(s_filtered_matrix2[i,j])) s_filtered_matrix2[i,j] <- mean(s_filtered_matrix2[i,(j-3):(j-1)])
    }
  }
  
  
  sw<-sw+1
  # Multiple recordings may be made within the same hour
  write.csv(s_matrix,file=paste("subject_items_at_dhours/",as.character(s),'_events.csv',sep=""))
  write.csv(s_filtered_matrix2,file=paste("subject_items_at_dhours/",as.character(s),'_filtered_events.csv',sep=""))
  print(paste('Subject ',as.character(sj),' ',as.character(sw),' profiling finished. hours: ',as.character(ncol(s_filtered_matrix2))))
}

