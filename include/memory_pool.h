/*
 * Copyright (C) 2014 Beat Küng <beat-kueng@gmx.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */

#ifndef _MEMORY_POOL_H_
#define _MEMORY_POOL_H_

#include <list>

/**
 ** class MemoryPool
 * memory pool based on linked lists. memory is only freed when this object is
 * deleted.
 */
template<class T>
class MemoryPool {
public:
	/**
	 * @param min_chunk_size   minimum number of free elements returned by next()
	 * @param alloc_factor     allocation chunk size = min_chunk_size * alloc_factor
	 */
	MemoryPool(int min_chunk_size, int alloc_factor = 100);
	~MemoryPool();

	/** reset: assume all memory is released */
	void reset();

	/** get next chunk (list of free elements), with at least min_chunk_size free elements */
	T* next();

	/** call this after next() to indicate how many elements were used */
	inline void setNumUsedElements(int num_used_elements);

	/** statistics: number of allocated chunks */
	std::size_t numChunks() const { return m_chunks.size(); }
private:
	void cleanup();
	void addNewChunk();

	struct Element {
		int length;
		int next_free;
		T* chunk;
	};
	std::list<Element> m_chunks;
	typename std::list<Element>::iterator m_cur_chunk;

	int m_min_chunk_size;
	int m_alloc_factor;
};

template<class T>
inline MemoryPool<T>::MemoryPool(int min_chunk_size, int alloc_factor)
	: m_min_chunk_size(min_chunk_size), m_alloc_factor(alloc_factor) {

	addNewChunk();
	reset();
}

template<class T>
inline MemoryPool<T>::~MemoryPool() {
	cleanup();
}

template<class T>
inline void MemoryPool<T>::reset() {
	m_cur_chunk = m_chunks.begin();
	m_cur_chunk->next_free = 0;
}
template<class T>
inline void MemoryPool<T>::addNewChunk() {
	Element new_element;
	new_element.length = m_min_chunk_size * m_alloc_factor;
	new_element.next_free = 0;
	new_element.chunk = new T[new_element.length];
	m_chunks.push_back(new_element);
}

template<class T>
T* MemoryPool<T>::next() {
	Element& element = *m_cur_chunk;
	if(element.length-element.next_free < m_min_chunk_size) {
		++m_cur_chunk;
		if(m_cur_chunk != m_chunks.end()) {
			//reuse next existing chunk
			m_cur_chunk->next_free = 0;
			return m_cur_chunk->chunk;
		}
		//allocate new chunk
		addNewChunk();
		m_cur_chunk = std::prev(m_chunks.end());
		return m_cur_chunk->chunk;
	}
	return m_cur_chunk->chunk + m_cur_chunk->next_free;
}

template<class T>
inline void MemoryPool<T>::setNumUsedElements(int num_used_elements) {
	m_cur_chunk->next_free += num_used_elements;
}

template<class T>
inline void MemoryPool<T>::cleanup() {
	for(auto& element : m_chunks) {
		delete[] element.chunk;
	}
	m_chunks.clear();
	m_cur_chunk = m_chunks.end();
}

#endif /* _MEMORY_POOL_H_ */
