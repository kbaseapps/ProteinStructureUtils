
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
 * <p>Original spec-file type: StructureToPDBFileParams</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "input_ref",
    "destination_dir"
})
public class StructureToPDBFileParams {

    @JsonProperty("input_ref")
    private String inputRef;
    @JsonProperty("destination_dir")
    private String destinationDir;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("input_ref")
    public String getInputRef() {
        return inputRef;
    }

    @JsonProperty("input_ref")
    public void setInputRef(String inputRef) {
        this.inputRef = inputRef;
    }

    public StructureToPDBFileParams withInputRef(String inputRef) {
        this.inputRef = inputRef;
        return this;
    }

    @JsonProperty("destination_dir")
    public String getDestinationDir() {
        return destinationDir;
    }

    @JsonProperty("destination_dir")
    public void setDestinationDir(String destinationDir) {
        this.destinationDir = destinationDir;
    }

    public StructureToPDBFileParams withDestinationDir(String destinationDir) {
        this.destinationDir = destinationDir;
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
        return ((((((("StructureToPDBFileParams"+" [inputRef=")+ inputRef)+", destinationDir=")+ destinationDir)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
