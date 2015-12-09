// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2015 enGits GmbH                                         +
// +                                                                    +
// + enGitsDataLib is free software: you can redistribute it and/or     +
// + modify it under the terms of the GNU Lesser General Public License +
// + as published by the Free Software Foundation, either version 3 of  +
// + the License, or (at your option) any later version.                +
// +                                                                    +
// + enGitsDataLib is distributed in the hope that it will be useful,   +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of     +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      +
// + GNU Lesser General Public License for more details.                +
// +                                                                    +
// + You should have received a copy of the GNU Lesser General Public   +
// + License along with enGitsDataLib.                                  +
// + If not, see <http://www.gnu.org/licenses/>.                        +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "list.h"
#include <iostream>
#include <stdlib.h>

#include <QStringList>

namespace EDL_NAMESPACE
{

void List::initList(size_t a_max_num_entries, size_t a_delta_entries)
{
  if (initialized) {
    //cout << "fatal error: attempt to initialize a List twice" << endl;
    //exit (EXIT_FAILURE);
    reset(a_max_num_entries, a_delta_entries);
  }
  m_MaxNumEntries = a_max_num_entries;
  m_DeltaEntries = a_delta_entries;
  //ATTENTION: Changed, July/01/04
  if(m_MaxNumEntries < 1) {
      m_MaxNumEntries = 1;
  }
  if (m_DeltaEntries < 1) {
      m_DeltaEntries = 1;
  }
  //END
  m_Offset = NULL;
  if (m_Master == this) {

    m_Active = new bool[m_MaxNumEntries];
    for (size_t i_entry = 0; i_entry < m_MaxNumEntries; i_entry++) {
      m_Active[i_entry] = false;
    }
    m_BlockSize = new size_t[m_MaxNumEntries];
    updateBlockSize();
    m_AfterLastEntry = 0;
    m_FirstBlock = 0;

    // create pthread-mutex
    //pthread_mutex_init(&m_Mutex, NULL);
  } else {
    m_Active = NULL;
    m_BlockSize = NULL;
  }
  m_Client = NULL;
  m_NumClients = 0;
  initialized = true;
}

void List::reset(size_t mne, size_t a_delta_entries)
{
  if (m_DeltaEntries == 0) {
    m_DeltaEntries = mne/10;
  }
  if (m_Master == this) {
    if (mne > m_MaxNumEntries) {
      extendList(mne - m_MaxNumEntries);
    }
    setDelta(a_delta_entries);
    for (unsigned int j = 0; j < m_NumClients; j++) {
      m_Client[j]->customReset(mne,a_delta_entries);
    }
  } else {
    m_Master->reset(mne, a_delta_entries);
  }
}

List::List()
{
  m_Client = NULL;
  m_Master = this;
  initialized = false;
  m_NumClients = 0;
  m_StateCounter = 0;
  m_Offset = NULL;
  m_Active = NULL;
  m_BlockSize = NULL;
  m_AfterLastEntry = 0;
  m_LinkNamesRequired = false;
  m_LinkLocked = false;
}

List::List(size_t a_max_num_entries, size_t a_delta_entries)
{
  m_IsClean = true;
  initialized = false;
  m_Master = this;
  initList (a_max_num_entries, a_delta_entries);
  m_NumClients = 0;
  m_LinkNamesRequired = false;
  m_LinkLocked = false;
}

List::List(List *linked_list, string link_name)
{
  m_Client = NULL;
  initialized = false;
  m_Master = linked_list->m_Master;
  initList (m_Master->maxNumEntries(), m_Master->deltaEntries());
  m_Master->addClient(this, link_name);
  m_NumClients = 0;
  m_LinkNamesRequired = false;
}

List::List(const List &other)
{
  m_Client = NULL;
  m_Master = this;
  initialized = false;
  m_NumClients = 0;
  m_StateCounter = 0;
  m_Offset = NULL;
  m_Active = NULL;
  m_BlockSize = NULL;
  m_AfterLastEntry = 0;
  operator=(other);
  m_LinkNamesRequired = false;
  m_LinkLocked = false;
}

void List::operator=(const List &other)
{
  if (initialized) {
    cerr << "fatal error: operator=() can only be udes for uninitialized Lists." << endl;
    exit(EXIT_FAILURE);
  }
  if (other.m_Master == &other) m_Master = this;
  else m_Master = other.m_Master;
  if (other.initialized) {
    initList(other.m_MaxNumEntries, other.m_DeltaEntries);
    // Fake that we are extending the list.
    // This is done to make sure that necessary memory will be allocated in derived classes
    size_t mne = m_MaxNumEntries;
    m_MaxNumEntries = 0;
    extend(mne);
    m_MaxNumEntries = mne;
  }
  if (this != m_Master) {
    m_Master->addClient(this);
  } else if (initialized) {
    m_AfterLastEntry = other.m_AfterLastEntry;
    for (size_t i = 0; i < m_MaxNumEntries; ++i) m_Active[i] = other.m_Active[i];
    updateBlockSize();
  }
}

void List::linkNamesOn(bool check_for_existing)
{
  if (this != m_Master) {
    cout << "'linkNammesOn() cannot be called for a 'List' which is already linked to another 'List'!" << endl;
    exit(EXIT_FAILURE);
  }
  if (m_NumClients > 0 && check_for_existing) {
    cout << "'This 'List' already has some link clients!" << endl;
    exit(EXIT_FAILURE);
  }
  m_LinkNamesRequired = true;
}

void List::linkNamesOff()
{
  if (this != m_Master) {
    cout << "'linkNammesOff() cannot be called for a 'List' which is already linked to another 'List'!" << endl;
    exit(EXIT_FAILURE);
  }
  m_LinkNamesRequired = false;
}

void List::lockLinking()
{
  if (this != m_Master) {
    master()->lockLinking();
  } else {
    m_LinkLocked = true;
  }
}

void List::unlockLinking()
{
  if (this != m_Master) {
    master()->unlockLinking();
  } else {
    m_LinkLocked = false;
  }
}

void List::link (List *linked_list, string link_name)
{
  if (this != m_Master) {
    cout << "This 'List' is already linked to another 'List'!" << endl;
    exit(EXIT_FAILURE);
  }
  if (linked_list->m_LinkLocked) {
    cout << "The 'List' has been locked for linking!" << endl;
    exit(EXIT_FAILURE);
  }
  if (m_NumClients != 0) {
    cout << "This 'List' is already a master for someone!" << endl;
    exit(EXIT_FAILURE);
  }
  m_Client = NULL;
  m_Master = linked_list->m_Master;
  if (master()->m_LinkNamesRequired) {
    if (link_name == "__none") {
      cout << "This 'List' requires link names!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  delAll();
  initList (m_Master->maxNumEntries(), m_Master->deltaEntries());

  // Fake that we are extending the list.
  // This is done to make sure that necessary memory will be allocated in derived classes
  size_t mne = m_MaxNumEntries;
  m_MaxNumEntries = 0;
  extend(mne);
  m_MaxNumEntries = mne;

  FORALL(i, this->) initEntry(i);
  m_Master->addClient(this, link_name);
}

void List::addClient (List *client_to_add, string link_name)
{
  if (this != m_Master) {
    cout << "fatal error trying to link to a non master list!" << endl;
    exit(EXIT_FAILURE);
  }
  if (m_LinkNamesRequired && link_name == "__none") {
    cout << "This 'List' requires link names!" << endl;
    exit(EXIT_FAILURE);
  }
  if (m_LinkNamesRequired && link_name != "__tmp") {
    for (size_t i = 0; i < m_NumClients; ++i) {
      if (m_Master->m_Client[i]->linkName() == link_name) {
        cout << "The name \"" << link_name.c_str() << "\" exists already in the ensemble of linked lists!" << endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  List **new_client;
  new_client = new List* [m_NumClients + 1];
  for (unsigned int i = 0; i < m_NumClients; i++) new_client[i] = m_Client[i];
  new_client[m_NumClients] = client_to_add;
  m_NumClients++;
  delete [] m_Client;
  m_Client = new_client;
  client_to_add->m_LinkName = link_name;
}

void List::delClient (List *client_to_del) {
  if (this != m_Master) {
    cout << "fatal error trying to delete a client from a non master list!" << endl;
    exit(EXIT_FAILURE);
  }
  List **new_client;
  int off = 0;
  for (unsigned int i = 0; i < m_NumClients - off; i++) {
    if (m_Client[i] == client_to_del) off += 1;
    m_Client[i] = m_Client[i + off];
  }
  new_client = new List* [m_NumClients - off];
  for (unsigned int i = 0; i < m_NumClients- off; i++) new_client[i] = m_Client[i];
  m_NumClients -= off;
  delete [] m_Client;
  m_Client = new_client;
  notifyChange();
}

List::~List () 
{
  if (this != m_Master) {
    m_Master->delClient(this);
  } else {
    for (size_t i = 0; i < m_NumClients; i++) {
      m_Client[i]->m_Master = m_Client[i];
    };

    // delete the mutex
    //pthread_mutex_destroy(&m_Mutex);

    delete [] m_BlockSize;
  }
  delete [] m_Offset;
  delete [] m_Client;
  delete [] m_Active;
}

size_t List::numActiveEntries () const
{
  size_t num_actives = 0;
  for (size_t i_entry = beginIdx(); i_entry < endIdx(); i_entry = nextIdx(i_entry)) {
    num_actives++;
  }
  return num_actives;
}

void List::extendList (size_t delta)
{
  if (this != m_Master) {
    cout << "fatal error in ExtendList this != master." << endl;
    exit(EXIT_FAILURE);
  }
  //cout << "void List::ExtendList(" << delta << ")" << endl;
  bool *new_active;
  extend (delta);
  m_MaxNumEntries += delta;
  new_active = new bool[m_MaxNumEntries];
  for (size_t i = 0; i < m_MaxNumEntries; i++) {
    if (i < m_MaxNumEntries - delta) {
      new_active[i] = m_Active[i];
    } else {
      new_active[i] = false;
    }
  }
  delete [] m_Active;
  m_Active = new_active;
  for (unsigned int j = 0; j < m_NumClients; j++) {
    m_Client[j]->extend(delta);
    m_Client[j]->m_MaxNumEntries += delta;
  }
  delete [] m_BlockSize;
  m_BlockSize = new size_t[m_MaxNumEntries];
  updateBlockSize();
}

void List::delEntry (size_t i_entry)
{
  if (!isActive(i_entry)) throw InvalidIndex_error(i_entry);
  for (size_t i_client = 0; i_client < m_NumClients; i_client++) {
    m_Client[i_client]->deleteData(i_entry);
  }
  deleteData(i_entry);
  m_Master->m_Active[i_entry] = false;
  m_Master->m_BlockSize[i_entry] = 1;
  m_IsClean = false;
  notifyChange();
}

void List::delAll () {
  if (initialized) {
    for (size_t i_entry = 0; i_entry < m_MaxNumEntries; i_entry++) {
      if (isActive(i_entry)) delEntry (i_entry);
    }
    cleanUp();
    notifyChange();
  }
}

void List::createOffsetList () {
  if (m_Master == this) {
    m_MaxOffset = 0;
    DeleteOffsetList();
    m_Offset = new size_t [endIdx()];
    for (size_t i_entry = 0; i_entry < endIdx(); i_entry++) {
      m_Offset[i_entry] = m_MaxOffset;
      if (!isActive(i_entry)) m_MaxOffset++;
    }
  } else {
    m_Master->createOffsetList();
  }
}

size_t List::offset (size_t i_entry) {
  if (m_Master->m_Offset == NULL) {
    cerr << "error: no offset information available" << endl;
    cerr << "       at this point in the program\n" << endl;
    exit (EXIT_FAILURE);
  } else {
    return m_Master->m_Offset[i_entry];
  }
  return 0; // dummy
}

void List::copy(size_t src, size_t dest)
{
  if (isMaster()) {
    copyEntry(src, dest);
    for (size_t i = 0; i < m_NumClients; i++) {
      m_Client[i]->copyEntry(src, dest);
    }
  } else {
    master()->copy(src, dest);
  }
}

void List::cleanUp () {
  size_t i_entry;
  if (m_Master == this) {
    createOffsetList();
    for (i_entry = 0; i_entry < endIdx(); i_entry++) {
      if (offset(i_entry) != 0) {
        copy(i_entry, i_entry - offset(i_entry));
      }
      if (i_entry < endIdx() - m_MaxOffset) {
        m_Active[i_entry] = true;
      } else {
        m_Active[i_entry] = false;
      }
    }
    updateBlockSize();
    try {
      //after_last_entry = prev_idx(max_num_entries) + 1;
      size_t i = m_MaxNumEntries;
      if (i == 0) {
        throw NotFound_error("List::prev");
      }
      do {
        i--;
      } while ((i > 0) && !isActive(i));
      if (!isActive(i)) {
        throw NotFound_error("List::prev");
      }
      m_AfterLastEntry = i + 1;
    } catch (NotFound_error) {
      //DBG_CP;
      m_AfterLastEntry = 0;
      //DBG_CP;
    }
    m_IsClean = true;
  } else {
    m_Master->cleanUp();
  }
}

void List::setDelta(size_t new_delta)
{
  if (this != m_Master) {
    cerr << "cannot modify the delta value for client lists" << endl;
    exit(EXIT_FAILURE);
  }
  m_DeltaEntries = new_delta;
}

void List::updateBlockSize()
{
  if (m_MaxNumEntries != 0) {
    m_BlockSize[m_MaxNumEntries - 1] = 1;
    for (size_t i = m_MaxNumEntries - 1; i > 0; i--) {
      if ((m_Active[i] && !m_Active[i-1]) || (!m_Active[i] && m_Active[i-1])) {
        m_BlockSize[i-1] = 1;
      } else {
        m_BlockSize[i-1] = m_BlockSize[i] + 1;
      }
    }
  }
  m_BlockSizeUpdated = true;
  m_FirstBlock = 0;
}

size_t List::tryAlloc(size_t n)
{
  if (isMaster()) {
    size_t new_entry = m_FirstBlock;
    bool found = false;
    do {
      if (!m_Active[new_entry] && (m_BlockSize[new_entry] >= n)) {
        found = true;
      } else {
        new_entry += m_BlockSize[new_entry];
      }
    } while ((new_entry < m_MaxNumEntries) && !found);
    if (found) {
      return new_entry;
    } else {
      if (m_BlockSizeUpdated) {
        extendList(m_DeltaEntries);
        m_BlockSizeUpdated = false;
      } else {
        updateBlockSize();
      }
      return tryAlloc(n);
    }
    return m_AfterLastEntry;
  } else {
    return m_Master->tryAlloc(n);
  }
}

size_t List::alloc(size_t n)
{
  if (isMaster()) {
    if (!initialized) {
      size_t d = n/10;
      if (d == 0) d = 1;
      initList(0,d);
      extendList(n);
    }
    size_t new_entry = tryAlloc(n);
    for (size_t i = new_entry; i < new_entry + n; i++) {
      m_BlockSize[i] = n - i + new_entry;
      m_Active[i] = true;
      initEntry(i);
      for (unsigned int j = 0; j < m_NumClients; j++) {
        m_Client[j]->initEntry(i);
      }
    }
    if (new_entry+1 > m_AfterLastEntry) m_AfterLastEntry = new_entry + n;
    m_FirstBlock = new_entry;
    notifyChange();
    return new_entry;
  } else {
    return m_Master->alloc(n);
  }
}

bool List::hasSameStructure(List *other_list)
{
  if (beginIdx() != other_list->beginIdx()) return false;
  if (endIdx() != other_list->endIdx()) return false;
  if (maxNumEntries() != other_list->maxNumEntries()) return false;
  for (size_t i = 0; i < maxNumEntries(); i++) {
    if (isActive(i) != other_list->isActive(i)) return false;
  }
  return true;
}

QStringList List::toBuffer(QByteArray& buffer)
{
  QStringList exclude_names;
  exclude_names << "__tmp";
  exclude_names << "__none";
  buffer.clear();
  QStringList names;
  if (isMaster()) {
    names << "__MASTER";
    for (size_t i = 0; i < endIdx(); ++i) {
      buffer += partialBuffer(i);
    }
    for (size_t i_client = 0; i_client < m_NumClients; ++i_client) {
      if (!exclude_names.contains(QString(m_Client[i_client]->linkName().c_str()))) {
        names << m_Client[i_client]->linkName().c_str();
        for (size_t i = 0; i < endIdx(); ++i) {
          buffer += m_Client[i_client]->partialBuffer(i);
        }
      }
    }
  } else {
    names = master()->toBuffer(buffer);
  }
  return names;
}

void List::fromBuffer(QByteArray buffer, QStringList names)
{
  if (isMaster()) {
    size_t idx = 0;
    foreach (QString name, names) {
      List* list = NULL;
      if (name == "__MASTER") {
        list = this;
      } else {
        for (size_t i_client = 0; i_client < m_NumClients; ++i_client) {
          if (m_Client[i_client]->linkName() == qPrintable(name)) {
            list = m_Client[i_client];
            break;
          }
        }
      }
      if (list) {
        for (size_t i = 0; i < endIdx(); ++i) {
          QByteArray partial_buffer = buffer.mid(idx, list->dataLength());
          list->fromPartialBuffer(i, partial_buffer);
          idx += list->dataLength();
        }
      }
    }
  } else {
    master()->fromBuffer(buffer, names);
  }
}

} // namespace





