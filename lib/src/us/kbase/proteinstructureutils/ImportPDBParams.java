
package us.kbase.proteinstructureutils;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: ImportPDBParams</p>
 * <pre>
 * Input of the import_matrix_from_excel function
 * input_shock_id: file shock id
 * input_file_path: absolute file path
 * input_staging_file_path: staging area file path
 * structure_name: structure object name
 * workspace_name: workspace name for object to be saved to
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "input_shock_id",
    "input_file_path",
    "input_staging_file_path",
    "structure_name",
    "description",
    "workspace_name"
})
public class ImportPDBParams {

    @JsonProperty("input_shock_id")
    private String inputShockId;
    @JsonProperty("input_file_path")
    private String inputFilePath;
    @JsonProperty("input_staging_file_path")
    private String inputStagingFilePath;
    @JsonProperty("structure_name")
    private String structureName;
    @JsonProperty("description")
    private String description;
    @JsonProperty("workspace_name")
    private String workspaceName;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("input_shock_id")
    public String getInputShockId() {
        return inputShockId;
    }

    @JsonProperty("input_shock_id")
    public void setInputShockId(String inputShockId) {
        this.inputShockId = inputShockId;
    }

    public ImportPDBParams withInputShockId(String inputShockId) {
        this.inputShockId = inputShockId;
        return this;
    }

    @JsonProperty("input_file_path")
    public String getInputFilePath() {
        return inputFilePath;
    }

    @JsonProperty("input_file_path")
    public void setInputFilePath(String inputFilePath) {
        this.inputFilePath = inputFilePath;
    }

    public ImportPDBParams withInputFilePath(String inputFilePath) {
        this.inputFilePath = inputFilePath;
        return this;
    }

    @JsonProperty("input_staging_file_path")
    public String getInputStagingFilePath() {
        return inputStagingFilePath;
    }

    @JsonProperty("input_staging_file_path")
    public void setInputStagingFilePath(String inputStagingFilePath) {
        this.inputStagingFilePath = inputStagingFilePath;
    }

    public ImportPDBParams withInputStagingFilePath(String inputStagingFilePath) {
        this.inputStagingFilePath = inputStagingFilePath;
        return this;
    }

    @JsonProperty("structure_name")
    public String getStructureName() {
        return structureName;
    }

    @JsonProperty("structure_name")
    public void setStructureName(String structureName) {
        this.structureName = structureName;
    }

    public ImportPDBParams withStructureName(String structureName) {
        this.structureName = structureName;
        return this;
    }

    @JsonProperty("description")
    public String getDescription() {
        return description;
    }

    @JsonProperty("description")
    public void setDescription(String description) {
        this.description = description;
    }

    public ImportPDBParams withDescription(String description) {
        this.description = description;
        return this;
    }

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public ImportPDBParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((("ImportPDBParams"+" [inputShockId=")+ inputShockId)+", inputFilePath=")+ inputFilePath)+", inputStagingFilePath=")+ inputStagingFilePath)+", structureName=")+ structureName)+", description=")+ description)+", workspaceName=")+ workspaceName)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
