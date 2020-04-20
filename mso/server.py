from flask import Flask, jsonify, make_response, request
import json
from cddd.inference import InferenceServer
from mso.optimizer import MPPSOOptimizerManualScoring

inferenceServer = InferenceServer(port_frontend=5520, use_running=True)
app = Flask(__name__)

@app.route('/init_swarm/', methods=['POST'])
def init_swarm():
    data = json.loads(request.data)
    optimizer = MPPSOOptimizerManualScoring.from_query(
        init_smiles=data["init_smiles"],
        num_part=data["num_part"],
        num_swarms=data["num_swarms"],
        inference_model=inferenceServer)
    output = [swarm.to_dict() for swarm in optimizer.swarms]
    return jsonify(output)

@app.route('/next_step/', methods=['POST'])
def next_step():
    data = json.loads(request.data)
    optimizer = MPPSOOptimizerManualScoring.from_swarm_dicts(
        swarm_dicts=data["swarm_dicts"],
        inference_model=inferenceServer,
        phi1=data["phi1"],
        phi2=data["phi2"],
        phi3=data["phi3"]
    )
    swarms = optimizer.run_one_iteration(data["fitness"])
    output = [swarm.to_dict() for swarm in swarms]
    return jsonify(output)


if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0', port=8890, ssl_context=('/gpfs01/home/ggwaq/.opensll/cert.pem', '/gpfs01/home/ggwaq/.opensll/key.pem'))
