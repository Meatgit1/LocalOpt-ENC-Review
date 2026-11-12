/*!
\file InTranceLibImpl.cpp
\version 1.7
\brief Implementation of shared library code to interface with gui.
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 05.02.2006 (AO)

\b Purpose:
 This file contains the implementation of a DLL - version of InTrance.
 It uses a multithreading library called ZThread to let InTrance run
 in a seperate thread so that the JAVA GUI won't get stuck while InTrance
 is calculating.


\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2001-2006 Bernd Dachwald, Andreas Ohndorf
*/

// InTrance
#include <configInTrance.h>
#include <InTranceConst.h>

#ifndef EXTRA_OUTPUT
 //#define EXTRA_OUTPUT 1
#endif


// STL
#include <iostream>
#include <string>
#include <sstream>


// possibly include JNI - headers from Java - GUI
// (Java Native Interface - JNI)
// and define directives for use of OS-dependent multithreading libraries
#ifdef JAVA_GUI
 #include <intrancegui_InTranceGUIFrame.h>

 // The following constants are used to adress an JAVA array of boolean
 // values that express the actual use of a float value data array.
 const int PLOT_SHOULD_UPDATE = 0;
 const int PLOT_IS_UPDATING   = 1;
 const int PLOT_IS_UPDATED    = 2;
 const int PLOT_IS_PAINTING   = 3;
#endif // JAVA_GUI


#ifdef COMMONCPP
// #include <cc++\thread.h>      // Common C++ library
#endif //COMMONCPP


#if defined (COMMONCPP)
 #include <TQueue.h>
#else
// #error InTranceLibImpl.cpp => no multithreading library
#endif

#if defined (WINDOWS) || defined (_WINDOWS)
 #include <windows.h>
#endif

// InTrance optimizer program
#include <CInTrance.h>


extern "C"
{
#ifdef __BORLANDC__
  #pragma argsused
#endif

// possible GUI
// The following code segement will be compiled just in case a
// JAVA graphical user interface (GUI) is used to present InTrance data to
// users. There are a few other code lines that are constrained to
// conditional compiling.

#ifdef JAVA_GUI

CInTrance* InTrance = 0;  // InTrance object
EAStatusRec  EA_status;   // data exchange structure
SimStatusRec SIM_status;  // simulation data exchange structure


// The used data puffer for caching data elements to transfer to GUI
// InTrance continuously generates data but the GUI only calls for data
// in certain time intervals. So that no data is get lost it will be cached
// in the following data buffer (see TQueue.h for definition).
#if defined COMMONCPP
 TransferDataItemBuffer* dataBuffer = 0;
#endif




/*
 * Class:     intrancegui_InTranceInterface
 * Method:    deleteDLLInTrance
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_deleteDLLInTrance
  (JNIEnv *env, jobject obj)
{
  // wait for InTrance to terminate
  while ((InTrance!=0) && (InTrance->isRunning() == true)) Sleep(250);
  SafeDelete(InTrance);
#if defined COMMONCPP
  SafeDelete(dataBuffer);
#endif
}





/*
 * Class:     InTranceInterface
 * Method:    initDLLInTrance
 * Signature: (Ljava/lang/String;)V
 */
JNIEXPORT void JNICALL
Java_intrancegui_InTranceGUIFrame_initDLLInTrance(JNIEnv *env, jobject obj, jstring file_name)
{
  // initialize our data buffer for transfers to JAVA GUI
#if defined COMMONCPP
  if (dataBuffer == 0) dataBuffer = new TransferDataItemBuffer;
#endif

  // converting file name expressed as JAVA string to C style string
  const char* converted_file_name = env->GetStringUTFChars(file_name,0);

  SafeDelete(InTrance);                   // delete possible old InTrance object
  InTrance = new CInTrance(converted_file_name);               // create new one

  // if creation of InTrance object has not been successful
  // throw an exception to JAVA caller
  if (InTrance == 0)
  {
    // get JAVA class id of java.lang.RuntimeException to indicate an error
    jclass newExcCls = env->FindClass("java/lang/RuntimeException");

    // there should be on because it's part of JAVA standard, but you never know
    if (newExcCls == 0) return; // return, if we cannot find exception class

    // throw exception
    env->ThrowNew(newExcCls, "error during allocation of InTrance object");
  }
  else
  {
    InTrance->setEAStatusObject(&EA_status);
  }
}


/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    isDLLInTranceRunning
 * Signature: ()Z
 */
JNIEXPORT jboolean JNICALL Java_intrancegui_InTranceGUIFrame_isDLLInTranceRunning
(JNIEnv *env, jobject obj) {

	if (InTrance != 0) {
	  if (InTrance->isRunning() == true) return jboolean(JNI_TRUE);
	}

	return jboolean(JNI_FALSE);
}


/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    newData
 * Signature: ()Z
 */
JNIEXPORT jboolean JNICALL Java_intrancegui_InTranceGUIFrame_newData
(JNIEnv *env, jobject obj)
{
	if (dataBuffer != 0)
	{
	  if (dataBuffer->empty() != true) return jboolean(JNI_TRUE);
	}
	return jboolean(JNI_FALSE);
}

//-----------------------------------------------------------------------------
void CopyTrajToJava(
        JNIEnv*           env
 ,      jobject           obj
 ,      jobject           transfer_item
 ,      jclass            cls
 ,const string            ids_n
 ,const string            ids_x
 ,const string            ids_y
 ,const string            ids_z
 ,      Trajectory        traj)
{
  if (traj != 0)                                      // is there a trajectory
  {
    if (traj->empty() != true) // and is it not empty
    {
      // IDs of Java s/c trajectory arrays
      jfieldID m_x_id(0),m_y_id(0),m_z_id(0),next_free_id(0);

      // jfloatArray objects of Java s/c trajectory arrays
      jfloatArray float_x_array(0),float_y_array(0),float_z_array(0);

      jfloat* x_array_start = 0; // pointer to first element of x - coordinate array
      jfloat* y_array_start = 0; // same for y - coordinate array
      jfloat* z_array_start = 0; // and z - coordinate array

      // get field ID (with JAVA field name "m_target_next_free")
      next_free_id = env->GetFieldID(cls, ids_n.c_str() ,"I");

      // Java Native Interface stuff to obtain pointers to array data structures
      // for x- coordinate array
      m_x_id = env->GetFieldID(cls, ids_x.c_str(), "[F");
      if (m_x_id == 0) Output(string("InTranceGUIFrame.getData (C++ impl.): can't find m_target_x!\n"));
      else 
      {
        float_x_array = reinterpret_cast<jfloatArray>(env->GetObjectField(transfer_item, m_x_id));
        x_array_start = env->GetFloatArrayElements(float_x_array,0);
      }
      // and y- coordinate array
      m_y_id = env->GetFieldID(cls, ids_y.c_str(), "[F");
      if (m_y_id == 0) Output(string("InTranceGUIFrame.getData (C++ impl.): can't find m_target_y!\n"));
      else
      {
        float_y_array = reinterpret_cast<jfloatArray>(env->GetObjectField(transfer_item, m_y_id));
        y_array_start = env->GetFloatArrayElements(float_y_array,0);
      }

      // and z- coordinate array
      m_z_id = env->GetFieldID(cls, ids_z.c_str(), "[F");
      if (m_z_id == 0) Output(string("InTranceGUIFrame.getData (C++ impl.): can't find m_target_z!\n"));
      else
      {
        float_z_array = reinterpret_cast<jfloatArray>(env->GetObjectField(transfer_item, m_z_id));
        z_array_start = env->GetFloatArrayElements(float_z_array,0);
      }

      if (x_array_start != 0 && y_array_start != 0 && z_array_start != 0) // are the desired arrays accessible
      {
        std::vector<MathLib::Vector>::iterator first = traj->begin();
        std::vector<MathLib::Vector>::iterator last  = traj->end();

        jint n_traj_points = jint(last-first); // calculate number of trajectory points
        jint array_size = env->GetArrayLength(float_x_array); //  get JAVA array size
        jint i=0;
        if (n_traj_points < array_size) // assure, that there is enough space in array
        {
          while (first != last)         // copy trajectory data
          {
            x_array_start[i] = jfloat((*first)(0));
            y_array_start[i] = jfloat((*first)(1));
            z_array_start[i] = jfloat((*first)(2));
            ++i;
            ++first;  // proceed to next point
          }
        }

        // set number of trajectory points
        env->SetIntField(transfer_item, next_free_id, n_traj_points);

        // release arrays so that JAVA can move them through garbage collecting
        env->ReleaseFloatArrayElements(float_x_array, x_array_start, 0);
        env->ReleaseFloatArrayElements(float_y_array, y_array_start, 0);
        env->ReleaseFloatArrayElements(float_z_array, z_array_start, 0);
      } // if (x_array_start!=0 && y_array_start !=0 && z_array_start!=0)
    } // if (traj->empty() != true)
  } // if (traj != 0)
} 


//-----------------------------------------------------------------------------
void copyDataArray2AnglePlot(
                           JNIEnv *env
                          ,jobject obj
                          ,jobject transfer_item
                          ,jclass cls
                          ,DataArray data
                          ,string& destination_plot_data_object
                          ,int dont_use_last_items = 0)
{

  const char* class_string = "Lintrancegui/PlotData;";

    if ( data != 0) // is there a trajectory
	  {
      if (data->empty() != true) // and is it not empty
		  {

        // first, we need access to PlotData class object of InTranceSolutionTransfer
        // it's name is "m_clock_angle_data"

        jfieldID plotDataId = 0;
        plotDataId = env->GetFieldID(              // obtain field id
          cls
         ,destination_plot_data_object.c_str()
         ,class_string);

        // the object field itself
        jobject plotDataObject = 0;
        plotDataObject = env->GetObjectField(transfer_item, plotDataId);

        // class id
        jclass plotDataClass = 0;
        plotDataClass = env->GetObjectClass(plotDataObject);

        // access variable of steering angle data array
	      jfieldID data_array_id(0);           // steering angle history array ID
	      jfloatArray data_array(0);	// steering angle history array
	      jfloat* data_array_start = 0;    // pointer to first array element

        // number of data items in steering angle array
	      jfieldID n_data_id(0);

        // access variables of boolean flag array
        jfieldID flags_array_id(0);   // boolean flags array ID
        jbooleanArray flagsArrayObject(0);  // boolean flags array
        jboolean* flagsArray = 0;           // pointer to first array element
        jint flagsArraySize(0);             // size of flags array

        // get ID of JAVA array of boolean flags indicating data array state
        flags_array_id = env->GetFieldID(plotDataClass,"m_flags","[Z");

        if (flags_array_id == 0)
         Output(string("InTranceGUIFrame.getData (C++ impl.): can't find \
          boolean array \"m_flags\"!\n"));
        else
        {
          // get JAVA object of boolean array
  	      flagsArrayObject = 
            reinterpret_cast<jbooleanArray>(env->GetObjectField(plotDataObject,
             flags_array_id));

          // in case of success determine number of boolean elements
          if (flagsArrayObject != 0) 
            flagsArraySize = env->GetArrayLength(flagsArrayObject);

          // there must be at least 5 of them; see also definition of
          // constants at the beginning of this file (eg PLOT_IS_UPDATING)
          if (flagsArraySize>=4)
          {
            // determine begin of array
            flagsArray = env->GetBooleanArrayElements(flagsArrayObject, 0);

            if (flagsArray[PLOT_SHOULD_UPDATE] == JNI_TRUE)
            {

	            // get field ID (with JAVA field name "m_n_data")
              n_data_id = env->GetFieldID(plotDataClass, "m_n_data","I");

              plotDataClass = env->GetObjectClass(plotDataObject);

	            // get ID of steering angle alpha array
	            data_array_id = env->GetFieldID(plotDataClass, "m_data", "[F");

              // in case of failure output warning message and do not continue
	            if (data_array_id == 0) Output(string("InTranceGUIFrame.getData (C++ impl.): can't find m_data!\n"));
	            else
              {
                // otherwise get object of array of float values
    	          data_array = 
                  reinterpret_cast<jfloatArray>(env->GetObjectField(plotDataObject,
                   data_array_id));

                // in case of success determine number of boolean elements
                jint data_array_size = 0;
                if (data_array != 0) 
                  data_array_size = env->GetArrayLength(data_array);

                CDataArray::iterator first = data->begin();
                CDataArray::iterator last  = data->end();
                jint n_hist_points = jint(last-first); // calculate number of trajectory points
                if (dont_use_last_items > 0)
                {
                  if (n_hist_points > dont_use_last_items)
                  {
                    last -= dont_use_last_items;
                    n_hist_points -= dont_use_last_items;
                  }
                }
                if (data_array_size >= n_hist_points) {
                  // obtain C - style pointer to beginning of array
	                data_array_start = env->GetFloatArrayElements(data_array,0);
        
  	              if (data_array_start != 0) // is the desired array accessible
                  {
                    // signal updating of data array by setting flag to true
                    flagsArray[PLOT_IS_UPDATING] = JNI_TRUE;
			              jint i=0;
			              if (n_hist_points < data_array_size) // assure, that there is enough space in array
                    {
                      while (first != last)         // copy angle data
                      {
              			    data_array_start[i] = jfloat((*first));
			                  ++i;
                        ++first;  // proceed to next point
                      }
                    }

			              // set number of trajectory points
		                env->SetIntField(plotDataObject, n_data_id, n_hist_points);
  
  		              // release arrays so that JAVA can move them through garbage collecting
	  	              env->ReleaseFloatArrayElements(data_array, data_array_start, 0);

                    // set update flag back to false because all data have been transfered
                    flagsArray[PLOT_IS_UPDATING] = JNI_FALSE;

                    // tell JAVA GUI, that new data is to be processed by setting flag
                    flagsArray[PLOT_IS_UPDATED] = JNI_TRUE;
 
                  } // if (data_array_start!=0)
                } // if (data_array_size >= n_hist_points) 
                else Output(string("getData (C++ impl.): not data array to hold data is too small!\n"));
              } //if (data_array_id == 0) 
              env->ReleaseBooleanArrayElements(flagsArrayObject, flagsArray, 0);
            } //if (flagsArray[SHOULD_UPDATE] == JNI_TRUE)
          } // if (flagsArraySize>=4)          
          else cout << "InTranceGUIFrame.getData (C++ impl.)=> bool array size < 4!\n";
        } // if (angle_alpha_flags_id == 0)
      } // if (item.m_angle_alpha_history->empty() != true)
	  } // if (item.m_angle_alpha_history != 0)

} // of copyDataArray2AnglePlot



/*
 * Class:     intrancegui_InTranceInterface
 * Method:    getData
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_getData
  (JNIEnv *env, jobject obj, jobject transfer_item)
{
  // get class id of transfer object
  jclass cls = env->GetObjectClass(transfer_item);

  // obtain variable id
  jfieldID fid = env->GetFieldID(cls, "m_reproduction_message", "Ljava/lang/String;");

  // get data from buffer
  TransferDataItem item = dataBuffer->get();

  // new data ?
  if (   item.m_trajectory != 0 
	  || item.m_target_trajectory  != 0
	  || item.m_reproduction_message.length() > 0)
  {

    jfieldID perf_id = env->GetFieldID(cls, "m_repro_per_sec", "F");
    if (perf_id != 0)
      env->SetFloatField(transfer_item, perf_id, jfloat(item.m_repro_per_sec));

    // copy s/c trajectory
    //-------------------------------------------------------------------------
    CopyTrajToJava(env, obj, transfer_item, cls
      , "m_next_free", "m_x", "m_y", "m_z", item.m_trajectory);

    // copy s/c trajectory (approx. steering angles)
    //-------------------------------------------------------------------------
    CopyTrajToJava(env, obj, transfer_item, cls
      , "m_approx_next_free", "m_approx_x", "m_approx_y", "m_approx_z"
      , item.m_approx_traj);

    // copy target trajectory
    //-------------------------------------------------------------------------
    CopyTrajToJava(env, obj, transfer_item, cls
      , "m_target_next_free", "m_target_x", "m_target_y", "m_target_z"
      , item.m_target_trajectory);

    // copy steering clock angle history
    //-------------------------------------------------------------------------
    copyDataArray2AnglePlot(
       env
      ,obj
      ,transfer_item
      ,cls
      ,item.m_angle_alpha_history
      ,string("m_clock_angle_data")
      );   

    // copy approximation of steering clock angle history
    //-------------------------------------------------------------------------
    copyDataArray2AnglePlot(
       env
      ,obj
      ,transfer_item
      ,cls
      ,item.m_angle_alpha_history_approx
      ,string("m_clock_angle_approx_data")
    );

    // copy steering cone angle history
    //-------------------------------------------------------------------------
    copyDataArray2AnglePlot(
       env
      ,obj
      ,transfer_item
      ,cls
      ,item.m_angle_beta_history
      ,string("m_cone_angle_data")
    );

    // copy approximation of steering cone angle history
    //-------------------------------------------------------------------------
    copyDataArray2AnglePlot(
       env
      ,obj
      ,transfer_item
      ,cls
      ,item.m_angle_beta_history_approx
      ,string("m_cone_angle_approx_data")
    );

    // copy thrust history
    //-------------------------------------------------------------------------
    copyDataArray2AnglePlot(
       env
      ,obj
      ,transfer_item
      ,cls
      ,item.m_thrust_history
      ,string("m_thrust_data")
    );

    // transfer reproduction status message
    //-------------------------------------------------------------------------
	  if (fid != 0) {
        env->SetObjectField(transfer_item, fid, 
          env->NewStringUTF(item.m_reproduction_message.c_str()));
	  }

    item.destroy(); // after transfering all the data, clear item

  }

  // set population size field
  fid = env->GetFieldID(cls, "m_ea_pop_size", "I");
  if (fid != 0) env->SetIntField(transfer_item, fid, jint(EA_status.population_size ));
  // set mutation probability field
  fid = env->GetFieldID(cls, "m_ea_mut_prob", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.mutation_probability));
  // set epoch field
  fid = env->GetFieldID(cls, "m_ea_epoch", "I");
  if (fid != 0) env->SetIntField(transfer_item, fid, jint(EA_status.epoch));
  // set # of bad epochs field
  fid = env->GetFieldID(cls, "m_ea_bad_epochs", "I");
  if (fid != 0) env->SetIntField(transfer_item, fid, jint(EA_status.bad_epochs));
  // set search space scan epoch field
  fid = env->GetFieldID(cls, "m_ea_sss_epoch", "I");
  if (fid != 0) env->SetIntField(transfer_item, fid, jint(EA_status.sss_epoch));
  // set # of sss epochs field
  fid = env->GetFieldID(cls, "m_ea_sss_epochs", "I");
  if (fid != 0) env->SetIntField(transfer_item, fid, jint(EA_status.n_sss_epochs));
  // set hypercube size field
  fid = env->GetFieldID(cls, "m_ea_hc_size", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.hc_size));
  // set hypercube size field
  fid = env->GetFieldID(cls, "m_ea_hc_limit", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.hc_upper_limit));
  // set hypercube shrinkinf factor
  fid = env->GetFieldID(cls, "m_ea_hc_sf", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.hc_sf));
  // set hypercube initial size
  fid = env->GetFieldID(cls, "m_ea_hc_initial_size", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.hc_init_size));
  // final boundary conditions met limit
  fid = env->GetFieldID(cls, "m_ea_fbc_met_limit", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.fbc_met_limit));
  // population convergence limit 1 (FBCs not met)
  fid = env->GetFieldID(cls, "m_pop_conv_limit1", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.pop_conv_limit1));
  // population convergence limit 2 (FBCs met)
  fid = env->GetFieldID(cls, "m_pop_conv_limit2", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.pop_conv_limit2));
  // EA convergence limit 1 (FBCs not met)
  fid = env->GetFieldID(cls, "m_ea_conv_limit1", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.ea_conv_limit1));
  // EA convergence limit 2 (FBCs met)
  fid = env->GetFieldID(cls, "m_ea_conv_limit2", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.ea_conv_limit2));
  // fitness of best individual of epoch
  fid = env->GetFieldID(cls, "m_ea_best_fitness_epoch", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.best_fitness_epoch));
  // fitness of best individual ever
  fid = env->GetFieldID(cls, "m_ea_best_fitness_ever", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.best_fitness_ever));
  // current best fitness of epoch
  fid = env->GetFieldID(cls, "m_ea_current_best_fitness", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.current_best_fitness));
  // rel. improvement of last best individual
  fid = env->GetFieldID(cls, "m_ea_improvement", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.improvement));
  // current average population fitness
  fid = env->GetFieldID(cls, "m_ea_avg_fitness", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.avg_fitness));
  // total # of reproductions
  fid = env->GetFieldID(cls, "m_ea_repro_ctr", "J");
  if (fid != 0) env->SetLongField(transfer_item, fid, jlong(EA_status.reproduction));
  // # of epoch reproductions
  fid = env->GetFieldID(cls, "m_ea_epoch_repro_ctr", "J");
  if (fid != 0) env->SetLongField(transfer_item, fid, jlong(EA_status.epoch_reproduction));
  // reproduction rate
  fid = env->GetFieldID(cls, "m_ea_repro_rate", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.reproduction_rate));
  // gen mutation probability
  fid = env->GetFieldID(cls, "m_ea_gen_mut_prob", "F");
  if (fid != 0) env->SetFloatField(transfer_item, fid, jfloat(EA_status.gen_mutation_probability));


  std::ostringstream os;
  for (unsigned int i=0;i<item.m_fitness_values.size();++i) 
  {
    os.clear();os.str("m_f");os.seekp(0,ios::end);
    os <<fixed<<setprecision(0)<<noshowpos<< (i+1);
    fid = env->GetFieldID(cls, os.str().c_str(), "Ljava/lang/String;");
    if (fid != 0) 
    {
      env->SetObjectField(transfer_item, fid, 
        env->NewStringUTF(item.m_fitness_values[i].c_str()));
	  }
  }

  // simulation data
  for (unsigned int i=0;i<SIM_status.data.size();++i) 
  {
    os.clear();os.str("m_sim");os.seekp(0,ios::end);
    os <<fixed<<setprecision(0)<<noshowpos<< (i+1);
    fid = env->GetFieldID(cls, os.str().c_str(), "Ljava/lang/String;");
    if (fid != 0) 
    {
      env->SetObjectField(transfer_item, fid, 
        env->NewStringUTF(SIM_status.data[i].c_str()));
	  }
  }

}




/*
 * Class:     InTranceInterface
 * Method:    startDLLInTrance
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_startDLLInTrance
  (JNIEnv *env, jobject obj)
{
#if defined (COMMONCPP)
  if (InTrance != 0)
  {
    InTrance->start(); // start InTrance thread
  } else 
  {
#if defined DEBUG
    cout << "InTranceGUIFrame.startDLLInTrance (C++ impl.) => InTrance object is null!\n";cout.flush();
#endif //DEBUG
  }
#else // COMMONCPP
  cout << "InTranceGUIFrame.startDLLInTrance (C++ impl.) => "
       << "No multithreading library! Need CommonC++!\n";
  cout.flush();
#endif //COMMONCPP
}




/*
 * Class:     InTranceInterface
 * Method:    stopDLLInTrance
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_stopDLLInTrance
  (JNIEnv *env, jobject obj)
{
  InTrance->stop();
  while (InTrance->isRunning()) Sleep(250);
  SafeDelete(InTrance);
}


/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    activateClockAngleHistory
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_activateClockAngleHistory
  (JNIEnv *, jobject)
{
  if (InTrance != 0) InTrance->activateClockAngleHistory();
}

/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    deactivateClockAngleHistory
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_deactivateClockAngleHistory
  (JNIEnv *, jobject)
{
  if (InTrance != 0) InTrance->deactivateClockAngleHistory();
}

/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    activateClockAngleHistoryApproximation
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_activateClockAngleHistoryApproximation
  (JNIEnv *, jobject)
{
  if (InTrance != 0) InTrance->activateClockAngleHistoryApproximation();
}

/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    deactivateClockAngleHistoryApproximation
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_deactivateClockAngleHistoryApproximation
  (JNIEnv *, jobject)
{
  if (InTrance != 0) InTrance->deactivateClockAngleHistoryApproximation();
}

/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    activateConeAngleHistory
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_activateConeAngleHistory
  (JNIEnv *, jobject)
{
  if (InTrance != 0) InTrance->activateConeAngleHistory();
}

/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    deactivateConeAngleHistory
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_deactivateConeAngleHistory
  (JNIEnv *, jobject)
{
  if (InTrance != 0) InTrance->deactivateConeAngleHistory();
}


/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    activateConeAngleHistoryApproximation
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_activateConeAngleHistoryApproximation
  (JNIEnv *, jobject)
{
  if (InTrance != 0) InTrance->activateConeAngleHistoryApproximation();
}

/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    deactivateConeAngleHistoryApproximation
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_deactivateConeAngleHistoryApproximation
  (JNIEnv *, jobject)
{
  if (InTrance != 0) InTrance->deactivateConeAngleHistoryApproximation();
}


/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    applyApproximationParameters
 * Signature: (II)V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_applyApproximationParameters
  (JNIEnv *, jobject, 
   jint clockSegSize,
   jint clockPolyDeg,
   jint coneSegSize,
   jint conePolyDeg)
{
  if (InTrance != 0) 
  {
    InTrance->setClockAngleApproximationParameters(
      clockSegSize,
      clockPolyDeg);
    InTrance->setConeAngleApproximationParameters(
      coneSegSize,
      conePolyDeg
    );
  }
}


/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    activateThrustHistory
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_activateThrustHistory
  (JNIEnv *, jobject)
{
  if (InTrance != 0) InTrance->activateThrustHistory();
}

/*
 * Class:     intrancegui_InTranceGUIFrame
 * Method:    deactivateThrustHistory
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_intrancegui_InTranceGUIFrame_deactivateThrustHistory
  (JNIEnv *, jobject)
{
  if (InTrance != 0) InTrance->deactivateThrustHistory();
}


#endif //JAVA_GUI

#ifdef WINDOWS
  // pregenerated DLL main function
  int WINAPI DllMain( HINSTANCE hInst,
    unsigned long reason,
    void * lpReserved )
  {
    return 1;
  }

  // pregenerated DLL entry point
  int WINAPI DllEntryPoint( HINSTANCE hInst, unsigned long reason, void * lpReserved )
  {
    return DllMain( hInst, reason, lpReserved );
  }
#endif //WINDOWS
}
