// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef LIST_H
#define LIST_H

#include "edlerror.h"

namespace EDL_NAMESPACE
{
class  List;
struct NotFound_error;
struct InvalidIndex_error;
}

#include <cstddef>
#include <cstdlib>
#include <typeinfo>
#include <mutex>

namespace EDL_NAMESPACE
{

/**
 * This exception will be thrown if an item cannot be found in a List,
 * or one of it's derived classes.
 */
struct NotFound_error : public EdlError
{
  int code;
  NotFound_error(std::string msg) : EdlError(msg) { code = 0; }
  NotFound_error(std::string msg, int a_code)  : EdlError(msg) { code = a_code; }
};

/// This exception will be thrown if 'mdebug' is enabled and an index is invalid.
struct InvalidIndex_error : public EdlError
{
  int i;
  InvalidIndex_error(size_t an_i)  : EdlError("InvalidIndex_error") { i = an_i; }
};

/**
 * This is a dynamic list. It provides only the necessary handling
 * routines and control data, but no payload. Other Lists can be linked to ensure
 * a consistent length of different Lists. They will expand and shrink simultaneously.
 * Please note that List and its
 * derived classes are not thread-safe in all cases.
 * It is, however, planned to make them thread-save.
 */
class List
{
public:

  typedef long unsigned int state_t;

private:

  /// the first block where new entries might be.
  size_t m_FirstBlock;

  /// one after the last active_entry.
  size_t m_AfterLastEntry;

  /// the dynamical increment.
  size_t m_DeltaEntries;

  /// the maximal number of entries.
  size_t m_MaxNumEntries;

  /// array to store the entry mapping after garbage collection.
  size_t *m_Offset;

  /// the maximal offset value found.
  size_t m_MaxOffset;

  /// indicates if a certain entry is active or not.
  bool *m_Active;

  /// number of client lists if this is a master.
  unsigned int m_NumClients;

  /// the master List, which is managing the dynamics.
  List *m_Master;

  /// a field to store pointers to the clients.
  List **m_Client;

  /// flag to determine if the list is clean (no holes).
  bool m_IsClean;

  /// field for the block sizes.
  size_t *m_BlockSize;

  /// update the block_size field to represent current state.
  void updateBlockSize();

  /// switch to determine if the block_size field is updated.
  bool m_BlockSizeUpdated;

  /// an internal counter representing list states @see State
  state_t m_StateCounter;

  /**
   * a mutex used if any maintenance function is called.
   * This will only happen in the master List.
   */
  //pthread_mutex_t m_Mutex;
  std::mutex m_Mutex;

  /// name of this list (for linked lists only)
  std::string m_LinkName;

  /// flag to determine if link names are required
  bool m_LinkNamesRequired;

  /// flag to disable linking to this list
  bool m_LinkLocked;


protected:

  /**
   * initialize the List.
   * @param mne the initial maximum number of entries
   * @param delta a delta entries to extend or shrink the list
   */
  void initList (size_t mne, size_t delta);

  /**
   * extend or shrink the list.
   * @param delta the number of new entries, negative values will shrink the list
   */
  void extendList(size_t delta);

  /**
   * add a new client if this is a master list.
   * @param client_to_add a pointer to the new client list
   */
  void addClient (List *client_to_add, std::string link_name = "__none");

  /**
   * delete a client from this master.
   * @param client_to_del a pointer to the client, which shall be deleted
   */
  void delClient (List *client_to_del);

  /// flag to determine if the list has been initialized.
  bool initialized;

  /**
   * initialization of an entry.
   * This method does nothing for the base-class and can be overwritten in child-classes.
   * @param i the index of the entry
   */
  virtual void initEntry (size_t i) {}

  /**
   * extend the list.
   * This method has to be overwritten in child-classes to ensure a proper
   * working of the cleaning mechanism.
   * @param delta the increment
   */
  virtual void extend (size_t delta) {}

  /**
   * copy an entry.
   * This method has to be overwritten in child-classes to ensure a proper
   * working of the cleaning mechanism.
   * @param src the index of the source entry
   * @param dest the index of the destination entry
   */
  virtual void copyEntry (size_t src, size_t dest) {}

  /**
   * notify a change of the list structure.
   * This method has to be called, whenever the list structure has been altered.
   * It results in an increased state counter.
   */
  void notifyChange() { ++m_StateCounter; }

  /// lock the mutex
  //void lockMutex() { pthread_mutex_lock(&(m_Master->m_Mutex)); }
  void lockMutex() { m_Master->m_Mutex.lock(); }

  /// release the mutex
  //void unlockMutex() { pthread_mutex_unlock(&(m_Master->m_Mutex)); }
  void unlockMutex() { m_Master->m_Mutex.unlock(); }

  virtual void customReset(real mne, real delta) {}

  virtual size_t dataLength() { return 0; }

  // QT needs to be removed (implement edl::ByteArry)
  // virtual QByteArray partialBuffer(size_t) { return QByteArray(); }
  // virtual void       fromPartialBuffer(size_t, QByteArray) {}

public:
  /**
   * get the master of this List.
   * @return a pointer to the master List
   */
  List* master() { return m_Master; }

  /**
   * reset the List.
   * @param mne the initial maximum number of entries
   * @param delta a delta entries to extend or shrink the list
   */
  virtual void reset (size_t mne, size_t delta = 0);

  /**
   * check if this is a master List.
   * @return true if it is a master List
   */
  bool isMaster() { return this == m_Master; }

  /**
   * get the maximal number of entries.
   * @ return the maximal number of entries
   */
  size_t maxNumEntries() const { return m_MaxNumEntries; }

  /**
   * get the current number of entries.
   * @return the current number of entries
   */
  size_t numEntries() const;
  size_t size() const { return numEntries(); }

  /**
   * get the dynamical increment.
   * @return the dynamical increment
   */
  size_t deltaEntries() const { return m_DeltaEntries; }

  /**
   * set the dynamical increment.
   * @param new_delta the new increment
   */
  void setDelta(size_t new_delta);

  /**
   * constructor.
   * Creates an uninitialized List.
   */
  List ();

  /**
   * constructor.
   * Creates a new List as master and initializes it.
   * @param mne initial maximum number of entries
   * @param delta initial increment of entries
   */
  List (size_t mne, size_t delta);

  /**
   * constructor.
   * Creates a new list which is linked to an existing master and
   * initializes it.
   * @param a_master the master List to link the new List to
   */
  List (List* a_master, std::string link_name = "__none");

  /// copy constructor.
  List(const List &other);

  /// destructor.
  virtual ~List ();

  /**
   * get the active flag for an entry.
   * @param i index of the entry
   * @return the active flag
   */
  bool isActive(size_t i) const {
    return m_Master->m_Active[i];
  }

  /**
   * get the number of active entries.
   * @return the number of active entries
   */
  size_t numActiveEntries () const;

  /**
   * add an entry at the end of the list.
   * @return the index of the new entry
   */
  size_t addEntry () { return alloc(1); }

  /**
   * allocate a number of new entries.
   * @param n the number of new entries
   * @return the index of the first new entry
   */
  size_t alloc(size_t n);

  /**
   * find out which would be the first entry
   * if a block would be allocated.
   * @param n the number of new entries
   * @return the index of the first new entry
   */
  size_t tryAlloc(size_t n);// { return num_entries; };

  /**
   * delete an entry
   * @param i the index of the entry to be deleted
   */
  void delEntry (size_t i);

  /**
   * delete a number of continuous entries.
   * @param i the first entry to delete
   * @param n the number of entries to be deleted
   */
  void deleteEntries(size_t i, size_t n)
  {
    while (n > 0) {
      delEntry(i+n-1);
      n--;
    }
  }

  /**
   * delete the data of an entry.
   * This method has to be overwritten if a derived class
   * allocates memory for new entries.
   * @param i the index of the entry to be deleted
   */
  virtual void deleteData(size_t i) { m_Master->m_BlockSize[i] += 0; } // avoid compiler warning :-(

  /// delete all entries.
  void delAll ();

  /**
   * find the next active entry.
   * This method is a bit slow right now, but that
   * will be improved in the future. So please stick to this and do not try
   * to write a fast workaround.
   * @param i the index of the entry whose following entry has to be found
   * @return the index of the following entry
   */
  size_t nextIdx(size_t i) const;

  /**
   * find the previous active entry.
   * This method is a bit slow right now, but that
   * will be improved in the future. So please stick to this and do not try
   * to write a fast workaround.
   * @param i the index of the entry whose following entry has to be found
   * @return the index of the following entry
   */
  size_t prevIdx(size_t i) const;

  /**
   * find the first active entry.
   * @return the index of the first active entry
   */
  size_t beginIdx() const;

  /**
   * find the last active enty.
   * @return the index of the last active entry
   */
  size_t lastIdx() const;

  /**
   * specifies an index after the last entry.
   * @return an index after the last one
   */
  size_t endIdx() const;

  /// create the offset-list for the cleaning process.
  void createOffsetList ();

  /// delete the offset list.
  void DeleteOffsetList () {
    delete [] m_Offset;
    m_Offset = NULL;
  }

  /**
   * get the offset of an entry.
   * @param i the index of the entry
   * @return the offset of the entry
   */
  size_t offset (size_t i);

  /**
   * find out if the list has deleted entries or if it is clean.
   * @return true if it is clean
   */
  bool isClean() const { return m_Master->m_IsClean; }

  /**
   * clean the list.
   * This method eliminates all holes from the list and resets num_entries.
   * After CleanUp num_entries will be identical with the number of active entries.
   */
  virtual void cleanUp ();

  /**
   * links this List to another List.
   * Please note that this can only be done for a 'standalone' List.
   * As soon as the List is linked to another list this method will
   * abort with an error message.
   * @param linked_list the List to link to
   */
  void link (List *linked_list, std::string link_name = "__none");

  /**
   * copy an entry.
   * @param src the index of the source entry
   * @param dest the index of the destination entry
   */
  virtual void copy (size_t src, size_t dest);

  /**
   * check if an index is valid.
   * This method check if an index is valid in this List.
   * If the index is invalid an InvalidIndex_error will be
   * thrown.
   * @param i the index to check
   */
  void checkIndex(size_t i)
  {
    if (!isActive(i)) {
      std::cerr << "index " << i << " is invalid" << std::endl;
      throw InvalidIndex_error(i);
    }
  }

  /**
   * check if another List has the same structure than this one.
   * @param other_list the List to compare
   * @return true if the structure is the same
   */
  bool hasSameStructure(List *other_list);

  /**
   * Get the state counter of the list.
   * The state counter can be used to determine if the structure has changed
   * (e.g. entries have been allocated or deleted).
   * Whenever the List changes the state counter will be increased.
   * @return the state counter
   */
  state_t state() { return m_Master->m_StateCounter; }

  /**
   * @brief Get the link name of this list.
   * @return the link name of this list.
   */
  std::string linkName() { return m_LinkName; }

  /**
   * @brief Enforce the use of link names for this list.
   * This method can only be called for master lists.
   * @param check_for_existing If set to true an error will be created if the list already has clients.
   */
  void linkNamesOn(bool check_for_existing = true);

  /**
   * @brief Do not enforce the use of link names for this list.
   * This method can only be called for master lists.
   */
  void linkNamesOff();

  /**
   * @brief diable linking
   */
  void lockLinking();

  /**
   * @brief re-enable linking
   */
  void unlockLinking();

  // QT needs to be removed (implement edl::ByteArry)
  // QStringList toBuffer(QByteArray& buffer);
  // void        fromBuffer(QByteArray buffer, QStringList names);

  void operator=(const List &other);
};



inline size_t List::endIdx() const
{
  return m_Master->m_AfterLastEntry;
}

inline size_t List::beginIdx() const
{
  size_t i = 0;
  while ((i < endIdx()) && !isActive(i)) ++i;
  return i;
}

inline size_t List::nextIdx(size_t i) const
{
  do { i++; } while ((i < endIdx()) && !isActive(i));
  return i;
}

inline size_t List::prevIdx(size_t i) const
{
  if (i == 0) throw NotFound_error("List::prev");
  do { i--; } while ((i > 0) && !isActive(i));
  if (!isActive(i)) {
    throw NotFound_error("List::prev");
  }
  return i;
}

inline size_t List::lastIdx() const
{
  if (endIdx() == 0) {
    throw NotFound_error("List::last");
  }
  return endIdx() - 1;
}

inline size_t List::numEntries() const
{
  return endIdx();
}


} // namsepace

/**
 * a loop over all active entries of a List.
 * This macro expands to for(size_t I = L begin(); I < L end(); I = L next(I)).
 * @param I the index variable
 * @param L the List, including the member-operator (for example "this->")
 */
#define FORALL(I,L) for (size_t I = L beginIdx(); I < L endIdx(); I = L nextIdx(I))

#endif
