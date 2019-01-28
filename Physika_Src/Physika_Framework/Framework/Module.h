#pragma once
#include "Physika_Core/Platform.h"
#include <memory>
#include <vector>
#include <cassert>
#include <iostream>
#include "Framework/Object.h"
#include "Framework/DeviceContext.h"

namespace Physika
{
class Node;

class Module : public Object
{
	DECLARE_CLASS(Module)
public:
	Module();
	Module(Node* node);

	virtual ~Module(void);

	virtual bool initialize() { return false; }

	virtual bool execute() { return false; }

	virtual bool updateStates() { return true; }

	int getInputSize() { return (int)m_ins.size(); }
	int getOutputSize() { return (int)m_outs.size(); }

	void setInputSize(int n) { m_ins.resize(n); }
	void setOutputSize(int n) { m_outs.resize(n); }

	void setInput(int id, std::string sematic)
	{
		assert(id < m_ins.size() && id >= 0);
		m_ins[id] = sematic;
	}

	void setOutput(int id, std::string sematic)
	{
		assert(id < m_outs.size() && id >= 0);
		m_outs[id] = sematic;
	}

	std::string getInput(int id)
	{
		assert(id < m_ins.size() && id >= 0);
		return m_ins[id];
	}

	std::string getOutput(int id)
	{
		assert(id < m_outs.size() && id >= 0);
		return m_outs[id];
	}

	void setNode(Node* node) {
		m_node = node;
	}

	std::shared_ptr<DeviceContext> getContext();

public:
	virtual bool insertToNode(Node* node) { return false; }
	virtual bool deleteFromNode(Node* node) { return false; }

private:
	std::vector<std::string> m_ins;
	std::vector<std::string> m_outs;

	Node* m_node;
};
}